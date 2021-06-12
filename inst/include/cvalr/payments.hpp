#pragma once

#include <iterator>
#include <numeric>
#include <vector>

#include "cvalr/derivative.hpp"

namespace cvalr {

namespace payments {

/**
 * Discounted payment leg
 *
 * @details
 * Wraps a container with coefficients \f$ c \f$ and provides an inner-product-like method:
 * \f[
 *   c^\top \pmatrix{1 \\ l_{t_1} \\ \vdots \\ l_{t_n} } .
 * \f]
 * The coefficients \f$ c \f$ represent an affine-linear function of (possibly) converted derivative
 * losses \f$ l_{t_1}, \ldots, l_{t_n} \f$.
 */
template <typename _Deriv>
class dl_functor {
 public:
  using deriv_type = _Deriv;

  template <typename _ForwardIt, typename _DerivIn>
  dl_functor(_DerivIn&& deriv, _ForwardIt c_first, _ForwardIt c_last)
      : deriv_{std::forward<_DerivIn>(deriv)}, coeffs_{c_first, c_last} {}

  template <typename _ContainerIn, typename _DerivIn>
  dl_functor(_DerivIn&& deriv, _ContainerIn&& coeffs)
      : deriv_{std::forward<_DerivIn>(deriv)}, coeffs_(std::forward<_ContainerIn>(coeffs)) {}

  auto coeffs() const { return coeffs_; }
  auto deriv() const { return deriv_; }

  template <typename _ForwardIt>
  double operator()(_ForwardIt l_first, _ForwardIt l_last) const {
    if (std::distance(l_first, l_last) + 1 != coeffs_.size()) {
      throw std::length_error("wrong iterator distance");
    }
    return coeffs_.front() +
           std::inner_product(
               ++coeffs_.cbegin(), coeffs_.cend(), l_first, 0., std::plus<double>{},
               [&deriv = deriv_](const auto x, const auto y) { return x * deriv(y); });
  }

 private:
  const deriv_type deriv_;
  const std::vector<double> coeffs_;
};

/**
 * Discounted default leg
 *
 * @details
 * For discounted-legs of the form
 * \f[
 *   \sum_{i=1}^n d_{t_i} [g(l_{t_i}) - g(l_{t_{i-1}})] .
 * \f]
 *
 * The coefficients are
 * \f[
 *   c_{i}
 *    = \begin{cases}
 *        0 & i = 0 , \\
 *        d_{t_i} - d_{t_{i+1}} & 0 < i < n , \\
 *        d_{t_n} & i = n .
 *      \end{cases}
 * \f]
 */
template <typename _Deriv>
class ddl_functor : public dl_functor<_Deriv> {
 public:
  using deriv_type = typename dl_functor<_Deriv>::deriv_type;
  template <typename _ForwardIt, typename _DerivIn>
  ddl_functor(_DerivIn&& deriv, _ForwardIt df_first, _ForwardIt df_last)
      : dl_functor<deriv_type>{std::forward<deriv_type>(deriv),
                               __create_coeffs(df_first, df_last)} {}

 private:
  template <typename _ForwardIt>
  std::vector<double> __create_coeffs(_ForwardIt df_first, _ForwardIt df_last) const {
    auto out = std::vector<double>{};
    out.reserve(std::distance(df_first, df_last) + 1);
    out.push_back(0.);
    std::copy(df_first, df_last, std::back_inserter(out));
    std::adjacent_difference(out.crbegin(), --out.crend(), out.rbegin());

    return out;
  }
};

/**
 * Discounted premium leg
 *
 * @details
 * For discounted-legs of the form
 * \f[
 *   p {[u - l]} + c \sum_{i=1}^{n} d_{t_i} \Delta t_{i-1} \frac{u - l - f l_{t_i} - f
 * l_{t_{i-1}}}{2} \f]
 */
template <typename _Deriv>
class dpl_functor : public dl_functor<_Deriv> {
 public:
  using deriv_type = typename dl_functor<_Deriv>::deriv_type;

  template <typename _ForwardIt1, typename _ForwardIt2, typename _DerivIn>
  dpl_functor(_DerivIn&& deriv, _ForwardIt1 t_first, _ForwardIt2 t_last, _ForwardIt2 df_first)
      : dl_functor<deriv_type>{deriv, __create_coeffs(t_first, t_last, df_first, deriv)} {}

 private:
  template <typename _ForwardIt1, typename _ForwardIt2, typename _DerivIn>
  static std::vector<double> __create_coeffs(_ForwardIt1 t_first, _ForwardIt1 t_last,
                                             _ForwardIt2 df_first, _DerivIn&& deriv) {
    auto ddt = std::vector<double>{};
    ddt.reserve(std::distance(t_first, t_last));
    std::adjacent_difference(t_first, t_last, std::back_inserter(ddt));
    std::transform(ddt.cbegin(), ddt.cend(), df_first, ddt.begin(), std::multiplies<double>{});

    auto out = std::vector<double>(ddt.size() + 1);
    out.front() = (deriv.upper - deriv.lower) *
                  (deriv.upfront + deriv.coupon * std::accumulate(ddt.cbegin(), ddt.cend(), 0.));
    std::adjacent_difference(ddt.crbegin(), ddt.crend(), out.rbegin(), std::plus<double>{});
    std::transform(++out.cbegin(), out.cend(), ++out.begin(), [&deriv](const auto val) {
      return -val * deriv.coupon * deriv.nominal_scale / 2.;
    });

    return out;
  }
};

template <typename _Deriv>
class dtl_functor : public dl_functor<_Deriv> {
 public:
  using deriv_type = typename dl_functor<_Deriv>::deriv_type;

  template <typename _ForwardIt1, typename _ForwardIt2, typename _DerivIn>
  dtl_functor(_DerivIn&& deriv, _ForwardIt1 t_first, _ForwardIt2 t_last, _ForwardIt2 df_first)
      : dl_functor<deriv_type>{deriv, __create_coeffs(t_first, t_last, df_first, deriv)} {}

 private:
  template <typename _ForwardIt1, typename _ForwardIt2, typename _DerivIn>
  std::vector<double> __create_coeffs(_ForwardIt1 t_first, _ForwardIt1 t_last, _ForwardIt2 df_first,
                                      _DerivIn&& deriv) const {
    const auto ddl_coeffs =
        ddl_functor<deriv_type>{deriv, df_first,
                                std::next(df_first, std::distance(t_first, t_last))}
            .coeffs();
    const auto dpl_coeffs = dpl_functor<deriv_type>{deriv, t_first, t_last, df_first}.coeffs();
    auto out = std::vector<double>{};
    out.reserve(ddl_coeffs.size());
    std::transform(ddl_coeffs.cbegin(), ddl_coeffs.cend(), dpl_coeffs.cbegin(),
                   std::back_inserter(out), std::minus<double>{});

    return out;
  }
};

}  // namespace payments

}  // namespace cvalr
