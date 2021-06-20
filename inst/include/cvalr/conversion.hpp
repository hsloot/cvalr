#pragma once

#include <algorithm>
#include <iterator>
#include <vector>

namespace cvalr {

namespace conversion {

template <typename _ForwardIt1, typename _ForwardIt2, typename _OutputIt>
void dt2adcp(_ForwardIt1 v_begin, _ForwardIt1 v_end, _ForwardIt2 t_begin, _ForwardIt2 t_end,
                     _OutputIt out_begin) {
  auto v = std::vector<double>{v_begin, v_end};
  std::sort(v.begin(), v.end());
  const auto d = static_cast<double>(v.size());
  for (auto [it_t, it_out, it_v] = std::make_tuple(t_begin, out_begin, v.cbegin()); it_t != t_end;
       ++it_t, ++it_out) {
    it_v = std::upper_bound(it_v, v.cend(), *it_t);
    *it_out = std::distance(v.cbegin(), it_v) / d;
  }
}

}  // namespace conversion

}  // namespace cvalr
