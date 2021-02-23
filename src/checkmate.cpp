#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// # nocov start

// [[Rcpp::export(rng=false)]]
bool is_exqmatrix(const NumericMatrix &x, const double tol) {
  if (x.nrow() != x.ncol())
    return false;
  const auto n = static_cast<std::size_t>(x.nrow());
  for (auto i = std::size_t{0}; i < n; ++i) {
    for (auto j = std::size_t{0}; j < n; ++j) {
      if (i > j && x(i, j) != 0.)
        return false;
      else if (i == j && x(i, j) > 0.)
        return false;
      else if (i < j && x(i, j) < 0.)
        return false;
    }
    const auto row = x(i, _);
    if (std::fabs(std::accumulate(row.cbegin(), row.cend(), 0.)) >= tol)
      return false;
  }
  return true;
}

// # nocov end
