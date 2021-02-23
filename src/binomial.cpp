#include <Rcpp.h>
#include <rmolib/math/binomial_coefficient.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
double multiply_binomial_coefficient(const double x, const std::size_t n,
                                     const std::size_t k) {
  return rmolib::math::multiply_binomial_coefficient(x, n, k);
}

// [[Rcpp::export]]
NumericVector v_multiply_binomial_coefficient(const NumericVector &x,
                                            const std::size_t n,
                                            const std::size_t k) {
  auto out = NumericVector(no_init(x.size()));
  std::transform(x.cbegin(), x.cend(), out.begin(),
                 [n = n, k = k](const auto v) {
                   return rmolib::math::multiply_binomial_coefficient(v, n, k);
                 });

  return out;
}
