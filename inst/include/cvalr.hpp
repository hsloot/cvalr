#pragma once

#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>
#include <vector>

namespace cvalr {

// ------------------------------  Pricing                            ------------------------------

/*! @brief Compute the discounted default leg
 *
 * @param l_begin Forward iterator for the begin-iterator of the *adcp* container
 * @param l_end Forward iterator for the end-iterator of the *adcp* container
 * @param df_begin Forward iterator for the begin-iterator of the *discount factors* container
 *   (assumed to have the same size as the *adcp* container)
 * @param l_map A unary operator mapping the *adcp* to the *derivative losses*
 */
template <typename _ForwardIt1, typename _ForwardIt2, typename _UnaryOp>
inline double ddl(_ForwardIt1 l_begin, _ForwardIt1 l_end, _ForwardIt2 df_begin, _UnaryOp l_map) {
  using std::distance;

  const auto n = static_cast<std::size_t>(distance(l_begin, l_end));

  auto dtl = std::vector<double>{};  // dtl: delta transformed loss
  dtl.reserve(n);
  std::adjacent_difference(
      l_begin, l_end, std::back_inserter(dtl),
      [&l_map](const auto val, const auto acc) { return l_map(val) - l_map(acc); });
  return std::inner_product(++dtl.cbegin(), dtl.cend(), ++df_begin, double{0});
}

/*! @brief Compute the discounted default leg of a *portfolio CDS*
 *
 * @param l_begin Forward iterator for the begin-iterator of the *adcp* container
 * @param l_end Forward iterator for the end-iterator of the *adcp* container
 * @param df_begin Forward iterator for the begin-iterator of the *discount factors* container
 *   (assumed to have the same size as the *adcp* container)
 * @param recovery_rate A number between zero and one for the *recovery rate* of the *portfolio CDS*
 */
template <typename _ForwardIt1, typename _ForwardIt2>
inline double ddl_pcds(_ForwardIt1 l_begin, _ForwardIt1 l_end, _ForwardIt2 df_begin,
                       const double recovery_rate) {
  const auto l_map = [recovery_rate](const auto val) { return val * (1. - recovery_rate); };
  return ddl(l_begin, l_end, df_begin, l_map);
}

/*! @brief Compute the discounted default leg of a *CDO tranche*
 *
 * @param l_begin Forward iterator for the begin-iterator of the *adcp* container
 * @param l_end Forward iterator for the end-iterator of the *adcp* container
 * @param df_begin Forward iterator for the begin-iterator of the *discount factors* container
 *   (assumed to have the same size as the *adcp* container)
 * @param recovery_rate A number between zero and one for the *recovery rate* of the *CDO tranche*
 * @param lower, upper Ordered numbers between zero and one for the
 *   *lower and upper attachment points* of the *CDO tranche*
 */
template <typename _ForwardIt1, typename _ForwardIt2>
inline double ddl_cdo(_ForwardIt1 l_begin, _ForwardIt1 l_end, _ForwardIt2 df_begin,
                      const double recovery_rate, const double lower, const double upper) {
  const auto l_map = [recovery_rate, lower, upper](auto val) {
    return std::min(std::max((1. - recovery_rate) * val - lower, 0.), upper - lower);
  };
  return ddl(l_begin, l_end, df_begin, l_map);
}

/*! @brief Compute the expected discounted default leg
 *
 * @param etl_begin Forward iterator for the begin-iterator of the *expected derivative loss*
 *   container
 * @param etl_end Forward iterator for the end-iterator of the *expected derivative loss* container
 * @param df_begin Forward iterator for the begin-iterator of the *discount factors* container
 *   (assumed to have the same size as the *adcp* container)
 */
template <typename _ForwardIt1, typename _ForwardIt2>
inline double eddl(_ForwardIt1 etl_begin, _ForwardIt1 etl_end, _ForwardIt2 df_begin) {
  return ddl(etl_begin, etl_end, df_begin, [](const auto val) { return val; });
}

/*! @brief Calculate the discounted premium leg (with spread 1)
 *
 * @param l_begin Forward iterator for the begin-iterator of the *adcp* container
 * @param l_end Forward iterator for the end-iterator of the *adcp* container
 * @param t_begin Forward iterator for the begin-iterator of *times* container
 *   (assumed to have the same size as the *adcp* container)
 * @param df_begin Forward iterator for the begin-iterator of the *discount factors* container
 *   (assumed to have the same size as the *adcp* container)
 * @param l_map A unary operator mapping the *adcp* to the *derivative losses*
 * @param n_map A unary operator mapping the *derivative losses* to the *(mid) nominal values*
 */
template <typename _ForwardIt1, typename _ForwardIt2, typename _UnaryOp1, typename _UnaryOp2>
inline double dpl1(_ForwardIt1 l_begin, _ForwardIt1 l_end, _ForwardIt2 t_begin,
                   _ForwardIt2 df_begin, _UnaryOp1 l_map, _UnaryOp2 n_map) {
  const auto n = static_cast<std::size_t>(std::distance(l_begin, l_end));

  auto nm = std::vector<double>{};  // nm: nominals (mid)
  nm.reserve(n);
  std::adjacent_difference(l_begin, l_end, std::back_inserter(nm),
                           [&l_map, &n_map](auto val, auto acc) {
                             val = n_map(l_map(val));
                             acc = n_map(l_map(acc));
                             return val + 0.5 * (acc - val);
                           });

  auto dt = std::vector<double>{};  // dt: delta times
  dt.reserve(n);
  std::adjacent_difference(t_begin, std::next(t_begin, n), std::back_inserter(dt));

  std::transform(nm.cbegin(), nm.cend(), dt.cbegin(), nm.begin(), std::multiplies<double>{});

  return std::inner_product(++nm.cbegin(), nm.cend(), ++df_begin, double{0});
}

/*! @brief Calculate the discounted premium leg (with spread 1) for a *portfolio CDS*
 *
 * @param l_begin Forward iterator for the begin-iterator of the *adcp* container
 * @param l_end Forward iterator for the end-iterator of the *adcp* container
 * @param t_begin Forward iterator for the begin-iterator of *times* container
 *   (assumed to have the same size as the *adcp* container)
 * @param df_begin Forward iterator for the begin-iterator of the *discount factors* container
 *   (assumed to have the same size as the *adcp* container)
 * @param recovery_rate A number between zero and one for the *recovery rate* of the *portfolio CDS*
 */
template <typename _ForwardIt1, typename _ForwardIt2>
inline double dpl1_pcds(_ForwardIt1 l_begin, _ForwardIt1 l_end, _ForwardIt2 t_begin,
                        _ForwardIt2 df_begin, const double recovery_rate) {
  const auto l_map = [recovery_rate](const auto val) { return val * (1. - recovery_rate); };
  const auto n_map = [recovery_rate](const auto val) { return 1. - val / (1. - recovery_rate); };

  return dpl1(l_begin, l_end, t_begin, df_begin, l_map, n_map);
}

/*! @brief Calculate the discounted premium leg (with spread 1) for a *CDO tranche*
 *
 * @param l_begin Forward iterator for the begin-iterator of the *adcp* container
 * @param l_end Forward iterator for the end-iterator of the *adcp* container
 * @param t_begin Forward iterator for the begin-iterator of *times* container
 *   (assumed to have the same size as the *adcp* container)
 * @param df_begin Forward iterator for the begin-iterator of the *discount factors* container
 *   (assumed to have the same size as the *adcp* container)
 * @param recovery_rate A number between zero and one for the *recovery rate* of the *CDO tranche*
 * @param lower, upper Ordered numbers between zero and one for the
 *   *lower and upper attachment points* of the *CDO tranche*
 */
template <typename _ForwardIt1, typename _ForwardIt2>
inline double dpl1_cdo(_ForwardIt1 l_begin, _ForwardIt1 l_end, _ForwardIt2 t_begin,
                       _ForwardIt2 df_begin, const double recovery_rate, const double lower,
                       const double upper) {
  const auto l_map = [recovery_rate, lower, upper](auto val) {
    return std::min(std::max((1. - recovery_rate) * val - lower, 0.), upper - lower);
  };
  const auto n_map = [lower, upper](const auto val) { return (upper - lower) - val; };

  return dpl1(l_begin, l_end, t_begin, df_begin, l_map, n_map);
}

/*! @brief Calculate the expected discounted premium leg (with spread 1)
 *
 * @param l_begin Forward iterator for the begin-iterator of the *expected derivative losses*
 *   container
 * @param l_end Forward iterator for the end-iterator of the *expected derivative losses* container
 * @param t_begin Forward iterator for the begin-iterator of *times* container
 *   (assumed to have the same size as the *expected derivative losses* container)
 * @param df_begin Forward iterator for the begin-iterator of the *discount factors* container
 *   (assumed to have the same size as the *expected derivative losses* container)
 * @param n_map A unary operator mapping the *expected derivative losses* to the
 *   *(mid) expected nominal values*
 */
template <typename _ForwardIt1, typename _ForwardIt2, typename _UnaryOp>
inline double edpl1(_ForwardIt1 l_begin, _ForwardIt1 l_end, _ForwardIt2 t_begin,
                    _ForwardIt2 df_begin, _UnaryOp n_map) {
  return dpl1(
      l_begin, l_end, t_begin, df_begin, [](const auto val) { return val; }, n_map);
}

/*! @brief Calculate the expected discounted premium leg (with spread 1) for a *portfolio CDS*
 *
 * @param l_begin Forward iterator for the begin-iterator of the *expected derivative losses*
 *   container
 * @param l_end Forward iterator for the end-iterator of the *expected derivative losses* container
 * @param t_begin Forward iterator for the begin-iterator of *times* container
 *   (assumed to have the same size as the *expected derivative losses* container)
 * @param df_begin Forward iterator for the begin-iterator of the *discount factors* container
 *   (assumed to have the same size as the *expected derivative losses* container)
 * @param recovery_rate A number between zero and one for the *recovery rate* of the *portfolio CDS*
 */
template <typename _ForwardIt1, typename _ForwardIt2>
inline double edpl1_pcds(_ForwardIt1 l_begin, _ForwardIt1 l_end, _ForwardIt2 t_begin,
                         _ForwardIt2 df_begin, const double recovery_rate) {
  const auto l_map = [](const auto val) { return val; };
  const auto n_map = [recovery_rate](const auto val) { return 1. - val / (1. - recovery_rate); };

  return dpl1(l_begin, l_end, t_begin, df_begin, l_map, n_map);
}

/*! @brief Calculate the expected discounted premium leg (with spread 1) for a *CDO tranche*
 *
 * @param l_begin Forward iterator for the begin-iterator of the *expected derivative losses*
 *   container
 * @param l_end Forward iterator for the end-iterator of the *expected derivative losses* container
 * @param t_begin Forward iterator for the begin-iterator of *times* container
 *   (assumed to have the same size as the *expected derivative losses* container)
 * @param df_begin Forward iterator for the begin-iterator of the *discount factors* container
 *   (assumed to have the same size as the *expected derivative losses* container)
 * @param recovery_rate A number between zero and one for the *recovery rate* of the *CDO tranche*
 * @param lower, upper Ordered numbers between zero and one for the
 *   *lower and upper attachment points* of the *CDO tranche*
 */
template <typename _ForwardIt1, typename _ForwardIt2>
inline double edpl1_cdo(_ForwardIt1 l_begin, _ForwardIt1 l_end, _ForwardIt2 t_begin,
                        _ForwardIt2 df_begin, const double lower, const double upper) {
  const auto l_map = [](const auto val) { return val; };
  const auto n_map = [lower, upper](const auto val) { return (upper - lower) - val; };

  return dpl1(l_begin, l_end, t_begin, df_begin, l_map, n_map);
}

// ------------------------------  Conversion                         ------------------------------

/*! @brief Calculate the *average default counting process* value from a container of
 *   *default times*
 *
 * @param v_begin Forward iterator for the begin-iterator of the *default times* container
 * @param v_end Forward iterator for the end-iterator of the *default times* container
 * @param t Non-negative number for the timepoint
 */
template <typename _ForwardIt1>
inline double dt2adcp(_ForwardIt1 v_begin, _ForwardIt1 v_end, const double t) {
  auto avg_smaller_than_t = [t, i = 0.](const auto acc, const auto val) mutable {
    ++i;
    return ((i - 1.) * acc + static_cast<decltype(acc)>(val <= t)) / i;
  };
  return std::accumulate(v_begin, v_end, 0., avg_smaller_than_t);
}

/*! @brief Convert a container of *default times* into a container of
 *    *average default counting process* values
 *
 * @param v_begin Forward iterator for the begin-iterator of the *default times* container
 * @param v_end Forward iterator for the end-iterator of the *default times* container
 * @param t_begin Forward iterator for the begin-iterator of the *timepoints* container
 * @param t_end Forward iterator for the end-iterator of the *timepoints* container
 * @param out_begin Forward iterator for the begin-iterator of the *output* container
 *   (assumed to have the same size as the *times* container)
 */
template <typename _ForwardIt1, typename _ForwardIt2, typename _OutputIt>
void dt2adcp(_ForwardIt1 v_begin, _ForwardIt1 v_end, _ForwardIt2 t_begin, _ForwardIt2 t_end,
             _OutputIt out_begin) {
  std::transform(t_begin, t_end, out_begin,
                 [v_begin, v_end](const auto t) { return dt2adcp(v_begin, v_end, t); });
}

/*! @brief Calculate the *pricing equation present value* from a container of
 *    *average default counting process values*
 *
 * @param v_begin Forward iterator for the begin-iterator of the *default times* container
 * @param v_end Forward iterator for the end-iterator of the *default times* container
 * @param t_begin Forward iterator for the begin-iterator of the *timepoints* container
 * @param t_end Forward iterator for the end-iterator of the *timepoints* container
 * @param df_begin Forward iterator for the begin-iterator of the *discount factors* container
 *   (assumed to have the same size as the *times* container)
 * @param coupon Number for the coupon
 * @param upfront Number for the upfront payment
 * @param l_map A unary operator mapping the *adcp* to the *derivative losses*
 * @param n_map A unary operator mapping the *derivative losses* to the *(mid) nominal values*
 * @param u_map A unary operator mapping the *quoted upfront* to the *paid upfront*
 */
template <typename _ForwardIt1, typename _ForwardIt2, typename _UnaryOp1, typename _UnaryOp2,
          typename _UnaryOp3>
inline double adcp2peqpv(_ForwardIt1 v_begin, _ForwardIt1 v_end, _ForwardIt2 t_begin,
                         _ForwardIt2 df_begin, const double coupon, const double upfront,
                         _UnaryOp1 l_map, _UnaryOp2 n_map, _UnaryOp3 u_map) {
  const auto ddl = cvalr::ddl(v_begin, v_end, df_begin, l_map);
  const auto dpl1 = cvalr::dpl1(v_begin, v_end, t_begin, df_begin, l_map, n_map);

  return ddl - (u_map(upfront) + coupon * dpl1);
}

/*! @brief Calculate the *pricing equation present value* for a *portfolio CDS* from a container of
 *    *average default counting process values*
 *
 * @param v_begin Forward iterator for the begin-iterator of the *default times* container
 * @param v_end Forward iterator for the end-iterator of the *default times* container
 * @param t_begin Forward iterator for the begin-iterator of the *timepoints* container
 * @param t_end Forward iterator for the end-iterator of the *timepoints* container
 * @param df_begin Forward iterator for the begin-iterator of the *discount factors* container
 *   (assumed to have the same size as the *times* container)
 * @param coupon Number for the coupon
 * @param upfront Number for the upfront payment
 * @param recovery_rate A number between zero and one for the *recovery rate* of the *portfolio CDS*
 */
template <typename _ForwardIt1, typename _ForwardIt2>
inline double adcp2peqpv_pcds(_ForwardIt1 v_begin, _ForwardIt1 v_end, _ForwardIt2 t_begin,
                              _ForwardIt2 df_begin, const double coupon, const double upfront,
                              const double recovery_rate) {
  const auto l_map = [recovery_rate](const auto val) { return val * (1. - recovery_rate); };
  const auto n_map = [recovery_rate](const auto val) { return 1. - val / (1. - recovery_rate); };
  const auto u_map = [](const auto val) { return val; };

  return adcp2peqpv(v_begin, v_end, t_begin, df_begin, coupon, upfront, l_map, n_map, u_map);
}

/*! @brief Calculate the *pricing equation present value* for a *CDO tranche* from a container of
 *    *average default counting process values*
 *
 * @param v_begin Forward iterator for the begin-iterator of the *default times* container
 * @param v_end Forward iterator for the end-iterator of the *default times* container
 * @param t_begin Forward iterator for the begin-iterator of the *timepoints* container
 * @param t_end Forward iterator for the end-iterator of the *timepoints* container
 * @param df_begin Forward iterator for the begin-iterator of the *discount factors* container
 *   (assumed to have the same size as the *times* container)
 * @param coupon Number for the coupon
 * @param upfront Number for the upfront payment
 * @param recovery_rate A number between zero and one for the *recovery rate* of the *CDO tranche*
 * @param lower, upper Ordered numbers between zero and one for the
 *   *lower and upper attachment points* of the *CDO tranche*
 */
template <typename _ForwardIt1, typename _ForwardIt2>
inline double adcp2peqpv_cdo(_ForwardIt1 v_begin, _ForwardIt1 v_end, _ForwardIt2 t_begin,
                             _ForwardIt2 df_begin, const double coupon, const double upfront,
                             const double recovery_rate, const double lower, const double upper) {
  const auto l_map = [recovery_rate, lower, upper](auto val) {
    return std::min(std::max((1. - recovery_rate) * val - lower, 0.), upper - lower);
  };
  const auto n_map = [lower, upper](const auto val) { return (upper - lower) - val; };
  const auto u_map = [lower, upper](const auto val) { return (upper - lower) * val; };

  return adcp2peqpv(v_begin, v_end, t_begin, df_begin, coupon, upfront, l_map, n_map, u_map);
}

}  // namespace cvalr
