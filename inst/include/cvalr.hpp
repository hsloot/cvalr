#pragma once

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

namespace cvalr {

template <typename _ForwardIt1, typename _ForwardIt2>
inline double eddl(_ForwardIt1 l_start, _ForwardIt1 l_end,
                   _ForwardIt2 df_start) {
  using std::distance;

  const auto n = distance(l_start, l_end);
  auto ldiff = std::vector<double>{};
  ldiff.reserve(n);
  std::adjacent_difference(l_start, l_end, std::back_inserter(ldiff));
  return std::inner_product(++ldiff.cbegin(), ldiff.cend(), ++df_start,
                            double{0});
}

template <typename _ForwardIt1, typename _ForwardIt2, typename _ForwardIt3,
          typename _UnaryOp>
inline double edpl1(_ForwardIt1 l_start, _ForwardIt1 l_end,
                    _ForwardIt2 t_start, _ForwardIt3 df_start, _UnaryOp mapn) {
  using std::distance;
  using std::next;

  const auto n = distance(l_start, l_end);
  auto tdiff = std::vector<double>{};
  tdiff.reserve(n);
  std::adjacent_difference(t_start, next(t_start, n),
                           std::back_inserter(tdiff));
  auto enmid = std::vector<double>{};
  enmid.reserve(n);
  std::transform(l_start, l_end, std::back_inserter(enmid), mapn);
  std::adjacent_difference(
      enmid.cbegin(), enmid.cend(), enmid.begin(),
      [](const auto val, const auto acc) { return val + 0.5 * (acc - val); });
  std::transform(tdiff.cbegin(), tdiff.cend(), enmid.cbegin(), enmid.begin(),
                 std::multiplies<double>{});
  return std::inner_product(++enmid.cbegin(), enmid.cend(), ++df_start,
                            double{0});
}

} // namespace cvalr
