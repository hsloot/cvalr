
// clang-format off
#include <Rcpp.h>
#include <rmolib/random/r_engine.hpp> // must be included before <rmolib/*>
// clang-format on

#include "rmolib/algorithm/r_shuffle.hpp"
#include "rmolib/random/multivariate/cuadras_auge_distribution.hpp"
#include "rmolib/random/multivariate/markovian_exmo_distribution.hpp"
#include "rmolib/random/univariate/exponential_distribution.hpp"
#include "rmolib/random/univariate/discrete_distribution.hpp"
#include "rmolib/random/univariate/uniform_int_distribution.hpp"
#include "rmolib/random/univariate/uniform_real_distribution.hpp"
#include <cvalr.hpp>

using namespace Rcpp;

static const R_xlen_t C_CHECK_USR_INTERRUP = 100000;

using exponential_distribution = rmolib::random::exponential_distribution<double>;
using uniform_real_distribution = rmolib::random::uniform_real_distribution<double>;
using uniform_int_distribution = rmolib::random::uniform_int_distribution<std::size_t>;
using discrete_distribution =
    rmolib::random::discrete_distribution<std::size_t, double, uniform_real_distribution,
                                          uniform_int_distribution>;
using cuadras_auge_distribution =
    rmolib::random::cuadras_auge_distribution<double, exponential_distribution>;
using markovian_exmo_distribution =
    rmolib::random::markovian_exmo_distribution<double, exponential_distribution,
                                                uniform_int_distribution, discrete_distribution,
                                                rmolib::algorithm::shuffler>;

struct no_shuffler {
  template <typename _RandomAccessIterator, typename _Engine, typename _UniformIntDistribution>
  void operator()(_RandomAccessIterator first, _RandomAccessIterator last, _Engine &&engine,
                  _UniformIntDistribution &&dist) const {
    return;
  }
};

using markovian_exmo_adcp_distribution = rmolib::random::markovian_exmo_distribution<
    double, exponential_distribution, uniform_int_distribution, discrete_distribution, no_shuffler>;

// [[Rcpp::export]]
NumericMatrix Rcpp__rexmo_markovian_acdp(const std::size_t n, const NumericVector &times,
                                         const std::size_t d, const NumericVector &ex_intensities) {
  using dist_t = markovian_exmo_adcp_distribution;
  using parm_t = typename dist_t::param_type;

  auto engine = r_engine{};
  const auto parm = parm_t{d, ex_intensities.begin(), ex_intensities.end()};
  auto dist = dist_t{};

  auto m = static_cast<std::size_t>(times.size());
  auto out = NumericMatrix(no_init(n, m));
  for (auto k = R_xlen_t{0}; k < n; ++k) {
    if ((d * k) % C_CHECK_USR_INTERRUP == 0) Rcpp::checkUserInterrupt();
    const auto values = dist(engine, parm);
    cvalr::dt2adcp(values.cbegin(), values.cend(), times.cbegin(), times.cend(), out(k, _).begin());
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix Rcpp__rcamo_esm_adcp(const std::size_t n, const NumericVector &times,
                                   const std::size_t d, const double alpha, const double beta) {
  using dist_t = cuadras_auge_distribution;
  using parm_t = typename dist_t::param_type;

  auto engine = r_engine{};
  const auto parm = parm_t{d, alpha, beta};
  auto dist = dist_t{};

  auto m = static_cast<std::size_t>(times.size());
  auto out = NumericMatrix(no_init(n, m));
  for (auto k = R_xlen_t{0}; k < n; ++k) {
    if ((d * k) % C_CHECK_USR_INTERRUP == 0) Rcpp::checkUserInterrupt();
    const auto values = dist(engine, parm);
    cvalr::dt2adcp(values.cbegin(), values.cend(), times.cbegin(), times.cend(), out(k, _).begin());
  }

  return out;
}

bool is_s4(SEXP x) { return ::Rf_isS4(x); }

S4 as_s4(SEXP x) {
  if (!is_s4(x)) stop("Not an S4 object");
  return S4(x);
}

template <typename _MODistribution>
typename _MODistribution::param_type extract_param(SEXP x, const double scale);

template <>
typename cuadras_auge_distribution::param_type extract_param<cuadras_auge_distribution>(
    SEXP x, const double scale) {
  using parm_t = typename cuadras_auge_distribution::param_type;
  const auto model = as_s4(x);
  if (!model.is("CuadrasAugeExtMO2FParam")) stop("Not an CuadrasAugeExtMO2FParam");
  const auto dim = as<std::size_t>(model.slot("dim"));
  const auto lambda = as<double>(model.slot("lambda"));
  const auto nu = as<double>(model.slot("nu"));
  const auto alpha = lambda * (1. - nu) * scale;
  const auto beta = lambda * nu * scale;

  return parm_t{dim, alpha, beta};
}

template <>
typename markovian_exmo_distribution::param_type extract_param<markovian_exmo_distribution>(
    SEXP x, const double scale) {
  using parm_t = typename markovian_exmo_distribution::param_type;
  const auto model = as_s4(x);
  if (!model.is("ExMOParam")) stop("Not an ExMOParam");
  auto ex_intensities = as<std::vector<double>>(NumericVector(model.slot("ex_intensities")));
  std::transform(ex_intensities.cbegin(), ex_intensities.cend(), ex_intensities.begin(),
                 [scale](const auto val) { return scale * val; });

  return parm_t{static_cast<std::size_t>(ex_intensities.size()), ex_intensities.cbegin(),
                ex_intensities.cend()};
}

template <>
typename markovian_exmo_adcp_distribution::param_type extract_param<markovian_exmo_adcp_distribution>(
    SEXP x, const double scale) {
  using parm_t = typename markovian_exmo_adcp_distribution::param_type;
  const auto model = as_s4(x);
  if (!model.is("ExMOParam")) stop("Not an ExMOParam");
  auto ex_intensities = as<std::vector<double>>(NumericVector(model.slot("ex_intensities")));
  std::transform(ex_intensities.cbegin(), ex_intensities.cend(), ex_intensities.begin(),
                 [scale](const auto val) { return scale * val; });

  return parm_t{static_cast<std::size_t>(ex_intensities.size()), ex_intensities.cbegin(),
                ex_intensities.cend()};
}

template <typename _MODistribution>
typename _MODistribution::param_type extract_global_param(const List &models, const double scale) {
  if (models.cbegin() == models.cend()) stop("Empty list");
  return extract_param<_MODistribution>(models[0], scale);
}

template <typename _MODistribution>
std::vector<typename _MODistribution::param_type> extract_partition_params(const List &models,
                                                                           const double scale) {
  if (models.cbegin() == models.cend()) stop("Empty list");
  if (++(models.cbegin()) == models.cend()) stop("One element list");
  auto out = std::vector<typename _MODistribution::param_type>{};
  for (auto [it, out_it] = std::make_pair(++(models.cbegin()), std::back_inserter(out));
       it != models.cend(); ++it, ++out_it) {
    *out_it = extract_param<_MODistribution>(*it, scale);
  }

  return out;
}

template <typename _MODistribution>
NumericMatrix Rcpp__rh2exmo_dt(const std::size_t n, const double fraction, const List &models) {
  using dist_t = _MODistribution;

  auto engine = r_engine{};
  const auto global_parm = extract_global_param<dist_t>(models, fraction);
  const auto partition_parms = extract_partition_params<dist_t>(models, 1. - fraction);
  auto dist = dist_t{};

  const auto d = global_parm.dim();
  const auto partition_sum =
      std::accumulate(partition_parms.cbegin(), partition_parms.cend(), std::size_t{0},
                      [](const auto acc, const auto &val) { return acc + val.dim(); });
  if (!(d == partition_sum)) stop("Invalid partition");

  auto out = NumericMatrix(no_init(n, d));
  for (auto k = R_xlen_t{0}; k < n; ++k) {
    if ((d * k) % C_CHECK_USR_INTERRUP == 0) Rcpp::checkUserInterrupt();
    const auto global_values = dist(engine, global_parm);
    auto partition_values = std::vector<double>{};
    partition_values.reserve(d);
    auto it = std::back_inserter(partition_values);
    for (const auto &parm : partition_parms) {
      const auto values = dist(engine, parm);
      it = std::copy(values.cbegin(), values.cend(), it);
    }
    std::transform(global_values.cbegin(), global_values.cend(), partition_values.cbegin(),
                   out(k, _).begin(), [](const auto x, const auto y) { return std::min(x, y); });
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix Rcpp__rh2exmo_markovian_dt(const std::size_t n, const double fraction,
                                         const List &models) {
  return Rcpp__rh2exmo_dt<markovian_exmo_distribution>(n, fraction, models);
}

// [[Rcpp::export]]
NumericMatrix Rcpp__rh2excamo_esm_dt(const std::size_t n, const double fraction,
                                     const List &models) {
  return Rcpp__rh2exmo_dt<cuadras_auge_distribution>(n, fraction, models);
}

template <typename _MODistribution>
NumericMatrix Rcpp__rh2exmo_adcp(const std::size_t n, const NumericVector &times,
                                 const double fraction, const List &models) {
  using dist_t = _MODistribution;

  auto engine = r_engine{};
  const auto global_parm = extract_global_param<dist_t>(models, fraction);
  const auto partition_parms = extract_partition_params<dist_t>(models, 1. - fraction);
  auto dist = dist_t{};

  const auto d = global_parm.dim();
  const auto partition_sum =
      std::accumulate(partition_parms.cbegin(), partition_parms.cend(), std::size_t{0},
                      [](const auto acc, const auto &val) { return acc + val.dim(); });
  if (!(d == partition_sum)) stop("Invalid partition");

  auto m = static_cast<std::size_t>(times.size());
  auto out = NumericMatrix(no_init(n, m));
  for (auto k = R_xlen_t{0}; k < n; ++k) {
    if ((d * k) % C_CHECK_USR_INTERRUP == 0) Rcpp::checkUserInterrupt();
    const auto global_values = dist(engine, global_parm);
    auto partition_values = std::vector<double>{};
    partition_values.reserve(d);
    auto it = std::back_inserter(partition_values);
    for (const auto &parm : partition_parms) {
      const auto values = dist(engine, parm);
      it = std::copy(values.cbegin(), values.cend(), it);
    }
    auto values = std::vector<double>(d);
    std::transform(global_values.cbegin(), global_values.cend(), partition_values.cbegin(),
                   values.begin(), [](const auto x, const auto y) { return std::min(x, y); });

    cvalr::dt2adcp(values.cbegin(), values.cend(), times.cbegin(), times.cend(), out(k, _).begin());
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix Rcpp__rh2exmo_markovian_adcp(const std::size_t n, const NumericVector &times,
                                           const double fraction, const List &models) {
  return Rcpp__rh2exmo_adcp<markovian_exmo_adcp_distribution>(n, times, fraction, models);
}

// [[Rcpp::export]]
NumericMatrix Rcpp__rh2excamo_esm_adcp(const std::size_t n, const NumericVector &times,
                                       const double fraction, const List &models) {
  return Rcpp__rh2exmo_adcp<cuadras_auge_distribution>(n, times, fraction, models);
}
