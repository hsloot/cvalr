#include <Rcpp.h>

#include <cvalr.hpp>

using namespace Rcpp;

//' Portfolio CDS spread
//'
//' @param expected_losses The expected losses of the portfolio (accounting
//'   for recovery) at `times`
//' @param expected_nominals The expected nominal values of the portfolio
//'   (accounting for expected accrued interest) at `times`
//' @param times The times of the payment schedule
//' @param discount_factors The discount factors.
//'
//' @export
// [[Rcpp::export]]
double portfolio_cds_spread(const NumericVector& expected_losses,
                            const NumericVector& expected_nominals,
                            const NumericVector& times,
                            const NumericVector& discount_factors) {
  return cvalr::eddl(expected_losses.cbegin(), expected_losses.cend(),
                     discount_factors.cbegin()) /
         cvalr::edpl1(expected_nominals.cbegin(), expected_nominals.cend(),
                      times.cbegin(), discount_factors.cbegin());
}
