#include <Rcpp.h>

#include <cvalr.hpp>

using namespace Rcpp;

//' Portfolio CDS spread and CDO methods
//'
//' @param expected_losses The expected losses (not including recovered part)
//' @param times The times of the payment schedule
//' @param discount_factors The discount factors.
//' @param recovery_rate The recovery rate
//'
//' @name cdx
//' @export
// [[Rcpp::export]]
double portfolio_cds_coupon(const NumericVector &expected_losses,
                            const NumericVector &times,
                            const NumericVector &discount_factors,
                            const double recovery_rate) {
  const auto eddl =
      cvalr::eddl(expected_losses.cbegin(), expected_losses.cend(),
                  discount_factors.cbegin());
  const auto edpl1 =
      cvalr::edpl1(expected_losses.cbegin(), expected_losses.cend(),
                   times.cbegin(), discount_factors.cbegin(),
                   [recovery_rate = recovery_rate](const auto val) {
                     return 1 - val / (1 - recovery_rate);
                   });
  return eddl / edpl1;
}

//' @rdname cdx
//'
//' @param coupon The periodically paid coupon of the payment leg
//' @param upfront The upfront payment
//'
//' @export
// [[Rcpp::export]]
double portfolio_cds_upfront(const NumericVector &expected_losses,
                             const NumericVector &times,
                             const NumericVector &discount_factors,
                             const double recovery_rate, const double coupon) {
  const auto eddl =
      cvalr::eddl(expected_losses.cbegin(), expected_losses.cend(),
                  discount_factors.cbegin());
  const auto edpl1 =
      cvalr::edpl1(expected_losses.cbegin(), expected_losses.cend(),
                   times.cbegin(), discount_factors.cbegin(),
                   [recovery_rate = recovery_rate](const auto val) {
                     return 1 - val / (1 - recovery_rate);
                   });
  return eddl - (coupon * edpl1);
}

//' @rdname cdx
//'
//' @param upfront The upfront payment
//'
//' @export
// [[Rcpp::export]]
double portfolio_cds_equation(const NumericVector &expected_losses,
                              const NumericVector &times,
                              const NumericVector &discount_factors,
                              const double recovery_rate, const double coupon,
                              const double upfront) {
  const auto eddl =
      cvalr::eddl(expected_losses.cbegin(), expected_losses.cend(),
                  discount_factors.cbegin());
  const auto edpl1 =
      cvalr::edpl1(expected_losses.cbegin(), expected_losses.cend(),
                   times.cbegin(), discount_factors.cbegin(),
                   [recovery_rate = recovery_rate](const auto val) {
                     return 1 - val / (1 - recovery_rate);
                   });
  return eddl - (upfront + coupon * edpl1);
}

//' @rdname cdx
//'
//' @param lower Lower attachment of the CDO tranche
//' @param upper Upper attachment of the CDO tranche
//'
//' @export
// [[Rcpp::export]]
double cdo_upfront(const NumericVector &expected_losses,
                   const NumericVector &times,
                   const NumericVector &discount_factors, const double lower,
                   const double upper, const double coupon) {
  const auto eddl =
      cvalr::eddl(expected_losses.cbegin(), expected_losses.cend(),
                  discount_factors.cbegin());
  const auto edpl1 = cvalr::edpl1(
      expected_losses.cbegin(), expected_losses.cend(), times.cbegin(),
      discount_factors.cbegin(),
      [lower, upper](const auto val) { return upper - lower - val; });
  return (eddl - coupon * edpl1) / (upper - lower);
}

//' @rdname cdx
//'
//' @export
// [[Rcpp::export]]
double cdo_coupon(const NumericVector &expected_losses,
                  const NumericVector &times,
                  const NumericVector &discount_factors, const double lower,
                  const double upper) {
  const auto eddl =
      cvalr::eddl(expected_losses.cbegin(), expected_losses.cend(),
                  discount_factors.cbegin());
  const auto edpl1 = cvalr::edpl1(
      expected_losses.cbegin(), expected_losses.cend(), times.cbegin(),
      discount_factors.cbegin(),
      [lower, upper](const auto val) { return upper - lower - val; });
  return eddl / edpl1;
}

//' @rdname cdx
//'
//' @export
// [[Rcpp::export]]
double cdo_equation(const NumericVector &expected_losses,
                    const NumericVector &times,
                    const NumericVector &discount_factors, const double lower,
                    const double upper, const double coupon,
                    const double upfront) {
  const auto eddl =
      cvalr::eddl(expected_losses.cbegin(), expected_losses.cend(),
                  discount_factors.cbegin());
  const auto edpl1 = cvalr::edpl1(
      expected_losses.cbegin(), expected_losses.cend(), times.cbegin(),
      discount_factors.cbegin(),
      [lower, upper](const auto val) { return upper - lower - val; });
  return eddl - ((upper - lower) * upfront + coupon * edpl1);
}
