#include <Rcpp.h>
#include <cvalr.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
double Rcpp__portfolio_cds_coupon(const NumericVector &expected_losses, const NumericVector &times,
                                  const NumericVector &discount_factors,
                                  const double recovery_rate) {
  const auto eddl =
      cvalr::eddl(expected_losses.cbegin(), expected_losses.cend(), discount_factors.cbegin());
  const auto edpl1 = cvalr::edpl1_pcds(expected_losses.cbegin(), expected_losses.cend(),
                                       times.cbegin(), discount_factors.cbegin(), recovery_rate);

  return eddl / edpl1;
}

// [[Rcpp::export]]
double Rcpp__portfolio_cds_upfront(const NumericVector &expected_losses, const NumericVector &times,
                                   const NumericVector &discount_factors,
                                   const double recovery_rate, const double coupon) {
  const auto eddl =
      cvalr::eddl(expected_losses.cbegin(), expected_losses.cend(), discount_factors.cbegin());
  const auto edpl1 = cvalr::edpl1_pcds(expected_losses.cbegin(), expected_losses.cend(),
                                       times.cbegin(), discount_factors.cbegin(), recovery_rate);

  return eddl - (coupon * edpl1);
}

// [[Rcpp::export]]
double Rcpp__portfolio_cds_equation(const NumericVector &expected_losses,
                                    const NumericVector &times,
                                    const NumericVector &discount_factors,
                                    const double recovery_rate, const double coupon,
                                    const double upfront) {
  const auto eddl =
      cvalr::eddl(expected_losses.cbegin(), expected_losses.cend(), discount_factors.cbegin());
  const auto edpl1 = cvalr::edpl1_pcds(expected_losses.cbegin(), expected_losses.cend(),
                                       times.cbegin(), discount_factors.cbegin(), recovery_rate);

  return eddl - (upfront + coupon * edpl1);
}

// [[Rcpp::export]]
double Rcpp__cdo_upfront(const NumericVector &expected_losses, const NumericVector &times,
                         const NumericVector &discount_factors, const double lower,
                         const double upper, const double coupon) {
  const auto eddl =
      cvalr::eddl(expected_losses.cbegin(), expected_losses.cend(), discount_factors.cbegin());
  const auto edpl1 = cvalr::edpl1_cdo(expected_losses.cbegin(), expected_losses.cend(),
                                      times.cbegin(), discount_factors.cbegin(), lower, upper);

  return (eddl - coupon * edpl1) / (upper - lower);
}

// [[Rcpp::export]]
double Rcpp__cdo_coupon(const NumericVector &expected_losses, const NumericVector &times,
                        const NumericVector &discount_factors, const double lower,
                        const double upper) {
  const auto eddl =
      cvalr::eddl(expected_losses.cbegin(), expected_losses.cend(), discount_factors.cbegin());
  const auto edpl1 = cvalr::edpl1_cdo(expected_losses.cbegin(), expected_losses.cend(),
                                      times.cbegin(), discount_factors.cbegin(), lower, upper);

  return eddl / edpl1;
}

// [[Rcpp::export]]
double Rcpp__cdo_equation(const NumericVector &expected_losses, const NumericVector &times,
                          const NumericVector &discount_factors, const double lower,
                          const double upper, const double coupon, const double upfront) {
  const auto eddl =
      cvalr::eddl(expected_losses.cbegin(), expected_losses.cend(), discount_factors.cbegin());
  const auto edpl1 = cvalr::edpl1_cdo(expected_losses.cbegin(), expected_losses.cend(),
                                      times.cbegin(), discount_factors.cbegin(), lower, upper);

  return eddl - ((upper - lower) * upfront + coupon * edpl1);
}
