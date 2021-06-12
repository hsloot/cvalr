#pragma once

#include <algorithm>

namespace cvalr {

namespace derivative {

struct derivative {
  derivative(const double recovery_rate = 0., const double lower = 0., const double upper = 1.,
             const double coupon = 0., const double upfront = 0., const double nominal_scale = 1.)
      : recovery_rate{recovery_rate},
        lower{lower},
        upper{upper},
        coupon{coupon},
        upfront{upfront},
        nominal_scale{nominal_scale} {}

  const double recovery_rate;
  const double lower;
  const double upper;
  const double coupon;
  const double upfront;
  const double nominal_scale;
};

struct pcds_derivative : public derivative {
  pcds_derivative(const double recovery_rate = 0., const double coupon = 0.,
                  const double upfront = 0.)
      : derivative(recovery_rate, 0., 1., coupon, upfront, 1. / (1. - recovery_rate)) {}

  double operator()(const double x) const { return (1 - recovery_rate) * x; }
};

struct cdo_derivative : public derivative {
  cdo_derivative(const double recovery_rate = 0., const double lower = 0., const double upper = 1.,
                 const double coupon = 0., const double upfront = 0.)
      : derivative(recovery_rate, lower, upper, coupon, upfront, 1.) {}

  double operator()(const double x) const {
    return std::min(std::max((1 - recovery_rate) * x - lower, 0.), upper - lower);
  }
};

template <typename _Deriv>
struct expectation_facade : public _Deriv {
  using _Deriv::_Deriv;
  double operator()(const double x) const { return x; }
};

using expected_pcds_derivative = expectation_facade<pcds_derivative>;
using expected_cdo_derivative = expectation_facade<cdo_derivative>;

}  // namespace derivative

}  // namespace cvalr
