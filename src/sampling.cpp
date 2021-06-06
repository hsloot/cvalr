
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

// [[Rcpp::export]]
NumericMatrix Rcpp__rexmo_markovian_acdp(const std::size_t n, const NumericVector &times,
                                         const std::size_t d, const NumericVector &ex_intensities) {
  using exponential_distribution = rmolib::random::exponential_distribution<double>;
  using uniform_real_distribution = rmolib::random::uniform_real_distribution<double>;
  using uniform_int_distribution = rmolib::random::uniform_int_distribution<std::size_t>;
  using discrete_distribution =
      rmolib::random::discrete_distribution<std::size_t, double, uniform_real_distribution,
                                            uniform_int_distribution>;
  using markovian_exmo_distribution =
      rmolib::random::markovian_exmo_distribution<double, exponential_distribution,
                                                  uniform_int_distribution, discrete_distribution,
                                                  rmolib::algorithm::shuffler>;

  using dist_t = markovian_exmo_distribution;
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
  using exponential_distribution = rmolib::random::exponential_distribution<double>;
  using cuadras_auge_distribution =
      rmolib::random::cuadras_auge_distribution<double, exponential_distribution>;

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
