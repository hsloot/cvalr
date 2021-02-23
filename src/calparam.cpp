#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(rng=false)]]
NumericMatrix dt2adcp(const NumericMatrix &x, const NumericVector &times) {
  const auto n = static_cast<std::size_t>(x.nrow());
  const auto d = static_cast<std::size_t>(x.ncol());
  const auto l = static_cast<std::size_t>(times.size());

  auto out = NumericMatrix(no_init(n, l));
  for (auto i = std::size_t{0}; i < n; ++i) {
    const auto values = x(i, _);
    for (auto j = std::size_t{0}; j < l; ++j) {
      out(i, j) = std::accumulate(
                      values.cbegin(), values.cend(), std::size_t{0},
                      [t = times[j]](const auto acc, const auto val) {
                        return acc + static_cast<decltype(acc)>(val <= t);
                      }) /
                  static_cast<double>(d);
    }
  }

  return out;
}

// [[Rcpp::export(rng=false)]]
NumericMatrix adcp2epd(const NumericMatrix &x, const std::size_t d) {
  const auto n = static_cast<std::size_t>(x.nrow());
  const auto l = static_cast<std::size_t>(x.ncol());

  auto out = NumericMatrix(no_init(d + 1, l));
  for (auto i = std::size_t{0}; i < d + 1; ++i) {
    for (auto j = std::size_t{0}; j < l; ++j) {
      const auto values = x(_, j);
      out(i, j) = std::accumulate(
          values.cbegin(), values.cend(), std::size_t{0},
          [v = static_cast<double>(i) /
                      static_cast<double>(d)](const auto acc, const auto val) {
            return acc + static_cast<decltype(acc)>(val == v);
          }) / static_cast<double>(n);
    }
  }

  return out;
}
