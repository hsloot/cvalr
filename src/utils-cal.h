#include <cstddef>
#include <numeric>
#include <utility>

template <typename _ForwardIterator1, typename _ForwardIterator2,
          typename _ForwardIterator3>
void dt2adcp(_ForwardIterator1 values_first, _ForwardIterator1 values_end,
             _ForwardIterator2 times_first, _ForwardIterator2 times_end,
             _ForwardIterator3 out_first) {
  const auto d =
      static_cast<std::size_t>(std::distance(values_first, values_end));

  for (auto [it_time, it_out] = std::make_pair(times_first, out_first);
       it_time != times_end; ++it_time, ++it_out) {
    *it_out =
        std::accumulate(values_first, values_end, std::size_t{0},
                        [t = *it_time](const auto acc, const auto val) {
                          return acc + static_cast<decltype(acc)>(val <= t);
                        }) /
        static_cast<double>(d);
  }
}
