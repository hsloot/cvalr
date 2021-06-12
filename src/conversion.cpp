#include <Rcpp.h>
#include <cvalr.hpp>
#include <fastcvalr.hpp>
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
  auto conversion = std::vector<cvalr::pcds_dtl_functor>{};
  for (auto j = std::size_t{0}; j < recovery_rate.size(); ++j) {
    conversion.emplace_back(cvalr::pcds_dtl_functor{times.cbegin(), times.cend(),
                                                    discount_factors.cbegin(), recovery_rate[j],
                                                    coupon[j], upfront[j]});
  }
  auto out = NumericMatrix(no_init(x.nrow(), recovery_rate.size()));
  for (auto k = size_t{0}; k < x.nrow(); ++k) {
    // copying is required because of issues with ConstMatrixRow's iterators
    const auto values = std::vector<double>(x(k, _).cbegin(), x(k, _).cend());
    for (auto j = size_t{0}; j < recovery_rate.size(); ++j) {
      out(k, j) = conversion[j](values.cbegin(), values.cend());
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
  auto conversion = std::vector<cvalr::cdo_dtl_functor>{};
  for (auto j = std::size_t{0}; j < recovery_rate.size(); ++j) {
    conversion.emplace_back(cvalr::cdo_dtl_functor{times.cbegin(), times.cend(),
                                                   discount_factors.cbegin(), recovery_rate[j],
                                                   lower[j], upper[j], coupon[j], upfront[j]});
  }
  auto out = NumericMatrix(no_init(x.nrow(), recovery_rate.size()));
  for (auto k = std::size_t{0}; k < x.nrow(); ++k) {
    // copying is required because of issues with ConstMatrixRow's iterators
    const auto values = std::vector<double>(x(k, _).cbegin(), x(k, _).cend());
    for (auto j = std::size_t{0}; j < recovery_rate.size(); ++j) {
      out(k, j) = conversion[j](values.cbegin(), values.cend());
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
      out(i, j) =
          std::min(std::max((1. - recovery_rate[j]) * x[i] - lower[j], 0.), upper[j] - lower[j]);
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
    out[i] = cvalr::pcds_edtl_functor{times.cbegin(),   times.cend(), discount_factors.cbegin(),
                                      recovery_rate[i], coupon[i],    upfront[i]}(values.cbegin(),
                                                                                  values.cend());
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
    out[i] = cvalr::cdo_edtl_functor{times.cbegin(),   times.cend(), discount_factors.cbegin(),
                                     recovery_rate[i], lower[i],     upper[i],
                                     coupon[i],        upfront[i]}(values.cbegin(), values.cend());
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
