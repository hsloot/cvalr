#include <Rcpp.h>
#include <cvalr.hpp>
using namespace Rcpp;

// [[Rcpp::export(rng=false)]]
NumericMatrix Rcpp__dt2adcp(const NumericMatrix &x, const NumericVector &times) {
  auto out = NumericMatrix(no_init(x.nrow(), times.size()));
  for (auto k = std::size_t{0}; k < x.nrow(); ++k) {
    const auto values = x(k, _);
    cvalr::dt2adcp(values.cbegin(), values.cend(), times.cbegin(), times.cend(), out(k, _).begin());
  }

  return out;
}

// [[Rcpp::export(rng=false)]]
NumericMatrix Rcpp__adcp2peqpv_pcds(const NumericMatrix &x, const NumericVector &times,
                                    const NumericVector &discount_factors,
                                    const NumericVector &recovery_rate, const NumericVector &coupon,
                                    const NumericVector &upfront) {
  auto out = NumericMatrix(no_init(x.nrow(), recovery_rate.size()));
  for (auto k = size_t{0}; k < x.nrow(); ++k) {
    // copying is required because of issues with ConstMatrixRow's iterators
    const auto values = std::vector<double>(x(k, _).cbegin(), x(k, _).cend());
    for (auto j = size_t{0}; j < recovery_rate.size(); ++j) {
      out(k, j) = cvalr::adcp2peqpv_pcds(values.cbegin(), values.cend(), times.cbegin(),
                                         discount_factors.cbegin(), coupon[j], upfront[j],
                                         recovery_rate[j]);
    }
  }
  return out;
}

// [[Rcpp::export(rng=false)]]
NumericMatrix Rcpp__adcp2peqpv_cdo(const NumericMatrix &x, const NumericVector &times,
                                   const NumericVector &discount_factors,
                                   const NumericVector &recovery_rate, const NumericVector &lower,
                                   const NumericVector &upper, const NumericVector &coupon,
                                   const NumericVector &upfront) {
  auto out = NumericMatrix(no_init(x.nrow(), recovery_rate.size()));
  for (auto k = std::size_t{0}; k < x.nrow(); ++k) {
    // copying is required because of issues with ConstMatrixRow's iterators
    const auto values = std::vector<double>(x(k, _).cbegin(), x(k, _).cend());
    for (auto j = std::size_t{0}; j < recovery_rate.size(); ++j) {
      out(k, j) = cvalr::adcp2peqpv_cdo(values.cbegin(), values.cend(), times.cbegin(),
                                        discount_factors.cbegin(), coupon[j], upfront[j],
                                        recovery_rate[j], lower[j], upper[j]);
    }
  }

  return out;
}

// [[Rcpp::export(rng=false)]]
NumericMatrix Rcpp__trans_v_pcds(const NumericVector &x, const NumericVector &recovery_rate) {
  auto out = NumericMatrix(no_init(x.size(), recovery_rate.size()));
  for (auto i = std::size_t{0}; i < x.size(); ++i) {
    const auto val = x[i];
    for (auto j = std::size_t{0}; j < recovery_rate.size(); ++j) {
      out(i, j) = (1 - recovery_rate[j]) * val;
    }
  }

  return out;
}

// [[Rcpp::export(rng=false)]]
NumericMatrix Rcpp__trans_v_cdo(const NumericVector &x, const NumericVector &recovery_rate,
                                const NumericVector &lower, const NumericVector &upper) {
  auto out = NumericMatrix(no_init(x.size(), recovery_rate.size()));
  for (auto i = std::size_t{0}; i < x.size(); ++i) {
    for (auto j = std::size_t{0}; j < recovery_rate.size(); ++j) {
      const auto val = (1 - recovery_rate[j]) * x[i];
      out(i, j) =
          (val < lower[j]) ? 0. : ((val > upper[j]) ? (upper[j] - lower[j]) : (val - lower[j]));
    }
  }

  return out;
}

// [[Rcpp::export(rng=false)]]
NumericVector Rcpp__lagg_ev_pcds(const NumericMatrix &x, const NumericVector &times,
                                 const NumericVector &discount_factors,
                                 const NumericVector &recovery_rate, const NumericVector &coupon,
                                 const NumericVector &upfront) {
  auto out = NumericVector(no_init(recovery_rate.size()));
  for (auto i = std::size_t{0}; i < recovery_rate.size(); ++i) {
    // copying is required because of issues with ConstMatrixRow's iterators
    const auto values = std::vector<double>(x(_, i).cbegin(), x(_, i).cend());
    const auto l_map = [](const auto val) { return val; };
    const auto n_map = [recovery_rate = recovery_rate[i]](const auto val) {
      return 1. - val / (1. - recovery_rate);
    };
    const auto u_map = [](const auto val) { return val; };

    out[i] =
        cvalr::adcp2peqpv(values.cbegin(), values.cend(), times.cbegin(), discount_factors.cbegin(),
                          coupon[i], upfront[i], l_map, n_map, u_map);
  }

  return out;
}

// [[Rcpp::export(rng=false)]]
NumericVector Rcpp__lagg_ev_cdo(const NumericMatrix &x, const NumericVector &times,
                                const NumericVector &discount_factors,
                                const NumericVector &recovery_rate, const NumericVector &lower,
                                const NumericVector &upper, const NumericVector &coupon,
                                const NumericVector &upfront) {
  auto out = NumericVector(no_init(recovery_rate.size()));
  for (auto i = std::size_t{0}; i < recovery_rate.size(); ++i) {
    // copying is required because of issues with ConstMatrixRow's iterators
    const auto values = std::vector<double>(x(_, i).cbegin(), x(_, i).cend());
    const auto l_map = [](const auto val) { return val; };
    const auto n_map = [lower = lower[i], upper = upper[i]](const auto val) {
      return (upper - lower) - val;
    };
    const auto u_map = [lower = lower[i], upper = upper[i]](const auto val) {
      return (upper - lower) * val;
    };
    out[i] =
        cvalr::adcp2peqpv(values.cbegin(), values.cend(), times.cbegin(), discount_factors.cbegin(),
                          coupon[i], upfront[i], l_map, n_map, u_map);
  }

  return out;
}

// [[Rcpp::export(rng=false)]]
NumericMatrix Rcpp__adcp2epd(const NumericMatrix &x, const std::size_t d) {
  auto out = NumericMatrix(no_init(d + 1, x.ncol()));
  for (auto i = std::size_t{0}; i < d + 1; ++i) {
    for (auto j = std::size_t{0}; j < x.ncol(); ++j) {
      out(i, j) =
          std::accumulate(x(_, j).cbegin(), x(_, j).cend(), 0.,
                          [v = static_cast<double>(i) / static_cast<double>(d), k = 0.](
                              const auto acc, const auto val) mutable {
                            ++k;
                            return ((k - 1.) * acc + static_cast<decltype(acc)>(val == v)) / k;
                          });
    }
  }

  return out;
}
