#pragma once

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

namespace cvalr {

template <typename _ForwardIterator>
inline double eddl(_ForwardIterator l_start, _ForwardIterator l_end,
                   _ForwardIterator d_start) {
  using std::distance;

  const auto n = distance(l_start, l_end);
  auto ldiff = std::vector<double>{};
  ldiff.reserve(n);
  std::adjacent_difference(l_start, l_end, std::back_inserter(ldiff));
  return std::inner_product(++ldiff.cbegin(), ldiff.cend(), ++d_start,
                            double{0});
}

template <typename _ForwardIterator>
inline double edpl1(_ForwardIterator n_start, _ForwardIterator n_end,
                    _ForwardIterator t_start, _ForwardIterator d_start) {
  using std::distance;
  using std::next;

  const auto n = distance(n_start, n_end);
  auto tdiff = std::vector<double>{};
  tdiff.reserve(n);
  std::adjacent_difference(t_start, next(t_start, n),
                           std::back_inserter(tdiff));
  auto nmid = std::vector<double>{};
  nmid.reserve(n);
  std::adjacent_difference(
      n_start, n_end, std::back_inserter(nmid),
      [](const auto val, const auto acc) { return val + 0.5 * (acc - val); });
  std::transform(tdiff.cbegin(), tdiff.cend(), nmid.cbegin(), nmid.begin(),
                 std::multiplies<double>{});
  return std::inner_product(++nmid.cbegin(), nmid.cend(), ++d_start, double{0});
}

}  // namespace cvalr
