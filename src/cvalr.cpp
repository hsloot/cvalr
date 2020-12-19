#include <Rcpp.h>

#include <cvalr.hpp>

using namespace Rcpp;

//' Portfolio CDS spread and CDO methods
//'
//' @param expected_default_counts The expected default counts of the portfolio
//'   at `times`
//' @param times The times of the payment schedule
//' @param discount_factors The discount factors.
//' @param recovery_rate The recovery rate
//'
//' @name cdx
//' @export
// [[Rcpp::export]]
double portfolio_cds_spread(const NumericVector& expected_default_counts,
                            const NumericVector& times,
                            const NumericVector& discount_factors,
                            const double recovery_rate) {
  const auto eddl = cvalr::eddl(expected_default_counts.cbegin(),
                                expected_default_counts.cend(),
                                discount_factors.cbegin(), recovery_rate);
  const auto edpl1 = cvalr::edpl1(expected_default_counts.cbegin(),
                                  expected_default_counts.cend(),
                                  times.cbegin(), discount_factors.cbegin());
  return eddl / edpl1;
}
