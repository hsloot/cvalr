#include <Rcpp.h>

#include <rmolib/math/binomial_coefficient.hpp>

using namespace Rcpp;

// [[Rcpp::export(rng=false)]]
NumericMatrix ex_intensities2qmatrix(const NumericVector& ex_intensities) {
  const auto d = static_cast<std::size_t>(ex_intensities.size());
  auto out = NumericMatrix(d + 1, d + 1);
  for (auto i = std::size_t{0}; i <= d; ++i) {
    if (i != d) {
      for (auto j = std::size_t{i + 1}; j <= d; ++j) {
        auto tmp = 0.;
        for (auto k = std::size_t{0}; k <= i; ++k) {
          tmp += rmolib::math::multiply_binomial_coefficient(
              ex_intensities[k + j - i - 1], i, k);
        }
        tmp = rmolib::math::multiply_binomial_coefficient(
            tmp, static_cast<std::size_t>(d - i),
            static_cast<std::size_t>(j - i));
        out(i, j) = tmp;
      }
      NumericMatrix::Row values = out(i, _);
      out(i, i) = -std::accumulate(values.cbegin(), values.cend(), 0.);
    }
  }
  return out;
}
