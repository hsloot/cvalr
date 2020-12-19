#pragma once

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

namespace cvalr {

template <typename _ForwardIt1, typename _ForwardIt2, typename _UnaryOp>
inline double eddl(_ForwardIt1 edc_start, _ForwardIt1 edc_end,
                   _ForwardIt2 df_start, _UnaryOp mapl) {
  using std::distance;

  const auto n = distance(edc_start, edc_end);
  auto eldiff = std::vector<double>{};
  eldiff.reserve(n);
  std::transform(edc_start, edc_end, std::back_inserter(eldiff), mapl);
  std::adjacent_difference(eldiff.cbegin(), eldiff.cend(), eldiff.begin());
  return std::inner_product(++eldiff.cbegin(), eldiff.cend(), ++df_start,
                            double{0});
}

template <typename _ForwardIt1, typename _ForwardIt2>
inline double eddl(_ForwardIt1 edc_start, _ForwardIt1 edc_end,
                   _ForwardIt2 df_start, const double recovery_rate) {
  return eddl(edc_start, edc_end, df_start, [recovery_rate](const auto val) {
    return (1 - recovery_rate) * val;
  });
}

template <typename _ForwardIt1, typename _ForwardIt2, typename _ForwardIt3,
          typename _UnaryOp>
inline double edpl1(_ForwardIt1 edc_start, _ForwardIt1 edc_end,
                    _ForwardIt2 t_start, _ForwardIt3 df_start, _UnaryOp mapn) {
  using std::distance;
  using std::next;

  const auto n = distance(edc_start, edc_end);
  auto tdiff = std::vector<double>{};
  tdiff.reserve(n);
  std::adjacent_difference(t_start, next(t_start, n),
                           std::back_inserter(tdiff));
  auto enmid = std::vector<double>{};
  enmid.reserve(n);
  std::transform(edc_start, edc_end, std::back_inserter(enmid), mapn);
  std::adjacent_difference(
      enmid.cbegin(), enmid.cend(), enmid.begin(),
      [](const auto val, const auto acc) { return val + 0.5 * (acc - val); });
  std::transform(tdiff.cbegin(), tdiff.cend(), enmid.cbegin(), enmid.begin(),
                 std::multiplies<double>{});
  return std::inner_product(++enmid.cbegin(), enmid.cend(), ++df_start,
                            double{0});
}

template <typename _ForwardIt1, typename _ForwardIt2, typename _ForwardIt3>
inline double edpl1(_ForwardIt1 edc_start, _ForwardIt1 edc_end,
                    _ForwardIt2 t_start, _ForwardIt3 df_start) {
  return edpl1(edc_start, edc_end, t_start, df_start,
               [](const auto val) { return 1 - val; });
}

}  // namespace cvalr
