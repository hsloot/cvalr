#pragma once

#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <Rcpp.h>

namespace cvalr {

class dl_functor {
 public:
  template <typename _ForwardIt>
  dl_functor(_ForwardIt first, _ForwardIt last) : coeffs_{first, last} {}

  dl_functor(std::vector<double>&& coeffs) : coeffs_(std::forward<std::vector<double>>(coeffs)) {}

  auto coeffs() const { return coeffs_; }

  // protected:
  template <typename _ForwardIt, typename _UnaryOp>
  double operator()(_ForwardIt l_first, _ForwardIt l_last, const _UnaryOp l_map) const {
    if (std::distance(l_first, l_last) + 1 != coeffs_.size()) {
      throw std::length_error("wrong iterator distance");
    }
    return coeffs_.front() +
           std::inner_product(++coeffs_.cbegin(), coeffs_.cend(), l_first, 0., std::plus<double>{},
                              [l_map](const auto x, const auto y) { return x * l_map(y); });
  }

 private:
  const std::vector<double> coeffs_;
};

class ddl_functor : public dl_functor {
 public:
  template <typename _ForwardIt>
  ddl_functor(_ForwardIt df_first, _ForwardIt df_last)
      : dl_functor{__create_coeffs(df_first, df_last)} {}

 private:
  template <typename _ForwardIt>
  static std::vector<double> __create_coeffs(_ForwardIt df_first, _ForwardIt df_last) {
    // TODO: Check range
    auto out = std::vector<double>{};
    out.reserve(std::distance(df_first, df_last) + 1);
    out.push_back(0.);
    std::copy(df_first, df_last, std::back_inserter(out));
    std::adjacent_difference(out.crbegin(), --out.crend(), out.rbegin());

    return out;
  }
};

class pcds_ddl_functor : public ddl_functor {
 public:
  template <typename _ForwardIt>
  pcds_ddl_functor(_ForwardIt df_first, _ForwardIt df_last, const double recovery_rate)
      : ddl_functor(df_first, df_last), recovery_rate_{recovery_rate} {}

  template <typename _ForwardIt>
  double operator()(_ForwardIt l_first, _ForwardIt l_last) const {
    const auto l_map = [recovery_rate = recovery_rate_](const auto val) {
      return (1. - recovery_rate) * val;
    };

    return dl_functor::operator()(l_first, l_last, l_map);
  }

 private:
  const double recovery_rate_;
};

class cdo_ddl_functor : public ddl_functor {
 public:
  template <typename _ForwardIt>
  cdo_ddl_functor(_ForwardIt df_first, _ForwardIt df_last, const double recovery_rate,
                  const double lower, const double upper)
      : ddl_functor(df_first, df_last),
        recovery_rate_{recovery_rate},
        lower_{lower},
        upper_{upper} {}

  template <typename _ForwardIt>
  double operator()(_ForwardIt l_first, _ForwardIt l_last) const {
    const auto l_map = [recovery_rate = recovery_rate_, lower = lower_,
                        upper = upper_](const auto val) {
      return std::min(std::max((1 - recovery_rate) * val - lower, 0.), upper - lower);
    };

    return dl_functor::operator()(l_first, l_last, l_map);
  }

 private:
  double recovery_rate_;
  double lower_;
  double upper_;
};

class eddl_functor : public ddl_functor {
 public:
  template <typename _ForwardIt>
  eddl_functor(_ForwardIt df_first, _ForwardIt df_last) : ddl_functor(df_first, df_last) {}

  template <typename _ForwardIt>
  double operator()(_ForwardIt l_first, _ForwardIt l_last) const {
    const auto l_map = [](const auto val) { return val; };

    return dl_functor::operator()(l_first, l_last, l_map);
  }
};

class dpl_functor : public dl_functor {
 public:
  template <typename _ForwardIt1, typename _ForwardIt2>
  dpl_functor(_ForwardIt1 t_first, _ForwardIt2 t_last, _ForwardIt2 df_first, const double lower,
              const double upper, const double coupon, const double upfront,
              const double nominal_scale)
      : dl_functor{__create_coeffs(t_first, t_last, df_first, lower, upper, coupon, upfront,
                                   nominal_scale)} {}

 private:
  template <typename _ForwardIt1, typename _ForwardIt2>
  static std::vector<double> __create_coeffs(_ForwardIt1 t_first, _ForwardIt1 t_last,
                                             _ForwardIt2 df_first, const double lower,
                                             const double upper, const double coupon,
                                             const double upfront, const double nominal_scale) {
    auto dt = std::vector<double>{};
    dt.reserve(std::distance(t_first, t_last));
    std::adjacent_difference(t_first, t_last, std::back_inserter(dt));

    auto ddt = std::vector<double>{};
    ddt.reserve(std::distance(t_first, t_last));
    std::transform(dt.cbegin(), dt.cend(), df_first, std::back_inserter(ddt),
                   std::multiplies<double>{});

    auto out = std::vector<double>(ddt.size() + 1);
    out.front() =
        (upper - lower) * (upfront + coupon * std::accumulate(ddt.cbegin(), ddt.cend(), 0.));
    std::adjacent_difference(ddt.crbegin(), ddt.crend(), out.rbegin(), std::plus<double>{});
    std::transform(
        ++out.cbegin(), out.cend(), ++out.begin(),
        [coupon, nominal_scale](const auto val) { return -val * coupon * nominal_scale / 2.; });

    return out;
  }
};

class pcds_dpl_functor : public dpl_functor {
 public:
  template <typename _ForwardIt1, typename _ForwardIt2>
  pcds_dpl_functor(_ForwardIt1 t_first, _ForwardIt1 t_last, _ForwardIt2 df_first,
                   const double recovery_rate, const double coupon, const double upfront)
      : dpl_functor(t_first, t_last, df_first, 0., 1., coupon, upfront, 1. / (1. - recovery_rate)),
        recovery_rate_{recovery_rate} {}

  template <typename _ForwardIt>
  double operator()(_ForwardIt l_first, _ForwardIt l_last) const {
    const auto l_map = [recovery_rate = recovery_rate_](const auto val) {
      return (1 - recovery_rate) * val;
    };

    return dl_functor::operator()(l_first, l_last, l_map);
  }

 private:
  const double recovery_rate_;
};

class cdo_dpl_functor : public dpl_functor {
 public:
  template <typename _ForwardIt1, typename _ForwardIt2>
  cdo_dpl_functor(_ForwardIt1 t_first, _ForwardIt1 t_last, _ForwardIt2 df_first,
                  const double recovery_rate, const double lower, const double upper,
                  const double coupon, const double upfront)
      : dpl_functor(t_first, t_last, df_first, lower, upper, coupon, upfront, 1.),
        recovery_rate_{recovery_rate},
        lower_{lower},
        upper_{upper} {}

  template <typename _ForwardIt>
  double operator()(_ForwardIt l_first, _ForwardIt l_last) const {
    const auto l_map = [recovery_rate = recovery_rate_, lower = lower_,
                        upper = upper_](const auto val) {
      return std::min(std::max((1 - recovery_rate) * val - lower, 0.), upper - lower);
    };

    return dpl_functor::operator()(l_first, l_last, l_map);
  }

 private:
  const double recovery_rate_;
  const double lower_;
  const double upper_;
};

class edpl_functor : public dpl_functor {
 public:
  template <typename _ForwardIt1, typename _ForwardIt2>
  edpl_functor(_ForwardIt1 t_first, _ForwardIt1 t_last, _ForwardIt2 df_first, const double lower,
               const double upper, const double coupon, const double upfront,
               const double nominal_scale)
      : dpl_functor(t_first, t_last, df_first, lower, upper, coupon, upfront, nominal_scale) {}

  template <typename _ForwardIt>
  double operator()(_ForwardIt l_first, _ForwardIt l_last) const {
    const auto l_map = [](const auto val) { return val; };

    return dpl_functor::operator()(l_first, l_last, l_map);
  }
};

class pcds_edpl_functor : public edpl_functor {
 public:
  template <typename _ForwardIt1, typename _ForwardIt2>
  pcds_edpl_functor(_ForwardIt1 t_first, _ForwardIt1 t_last, _ForwardIt2 df_first,
                    const double recovery_rate, const double coupon, const double upfront)
      : edpl_functor(t_first, t_last, df_first, 0., 1., coupon, upfront,
                     1. / (1. - recovery_rate)) {}
};

class cdo_edpl_functor : public edpl_functor {
 public:
  template <typename _ForwardIt1, typename _ForwardIt2>
  cdo_edpl_functor(_ForwardIt1 t_first, _ForwardIt1 t_last, _ForwardIt2 df_first,
                   const double recovery_rate, const double lower, const double upper,
                   const double coupon, const double upfront)
      : edpl_functor(t_first, t_last, df_first, lower, upper, coupon, upfront, 1.) {}
};

class dtl_functor : public dl_functor {
 public:
  dtl_functor(std::vector<double>&& ddl_coeffs, std::vector<double>&& dpl_coeffs)
      : dl_functor{__create_coeffs(std::forward<std::vector<double>>(ddl_coeffs),
                                   std::forward<std::vector<double>>(dpl_coeffs))} {}

 private:
  std::vector<double> __create_coeffs(std::vector<double>&& ddl_coeffs,
                                      std::vector<double>&& dpl_coeffs) const {
    auto out = std::vector<double>(std::forward<std::vector<double>>(ddl_coeffs));
    std::transform(out.cbegin(), out.cend(), dpl_coeffs.cbegin(), out.begin(), std::minus<double>{});

    return out;
  }
};

class pcds_dtl_functor : public dtl_functor {
 public:
  template <typename _ForwardIt1, typename _ForwardIt2>
  pcds_dtl_functor(_ForwardIt1 t_first, _ForwardIt1 t_last, _ForwardIt2 df_first,
                   const double recovery_rate, const double coupon, const double upfront)
      : dtl_functor(
            pcds_ddl_functor{df_first, std::next(df_first, std::distance(t_first, t_last)),
                             recovery_rate}
                .coeffs(),
            pcds_dpl_functor{t_first, t_last, df_first, recovery_rate, coupon, upfront}.coeffs()),
        recovery_rate_{recovery_rate} {}

  template <typename _ForwardIt>
  double operator()(_ForwardIt l_first, _ForwardIt l_last) const {
    const auto l_map = [recovery_rate = recovery_rate_](const auto val) {
      return (1 - recovery_rate) * val;
    };

    return dl_functor::operator()(l_first, l_last, l_map);
  }

 private:
  const double recovery_rate_;
};

class cdo_dtl_functor : public dtl_functor {
 public:
  template <typename _ForwardIt1, typename _ForwardIt2>
  cdo_dtl_functor(_ForwardIt1 t_first, _ForwardIt1 t_last, _ForwardIt2 df_first,
                  const double recovery_rate, const double lower, const double upper,
                  const double coupon, const double upfront)
      : dtl_functor(
            cdo_ddl_functor{df_first, std::next(df_first, std::distance(t_first, t_last)),
                            recovery_rate, lower, upper}
                .coeffs(),
            cdo_dpl_functor{t_first, t_last, df_first, recovery_rate, lower, upper, coupon, upfront}
                .coeffs()),
        recovery_rate_{recovery_rate},
        lower_{lower},
        upper_{upper} {}

  template <typename _ForwardIt>
  double operator()(_ForwardIt l_first, _ForwardIt l_last) const {
    const auto l_map = [recovery_rate = recovery_rate_, lower = lower_,
                        upper = upper_](const auto val) {
      return std::min(std::max((1 - recovery_rate) * val - lower, 0.), upper - lower);
    };

    return dl_functor::operator()(l_first, l_last, l_map);
  }

 private:
  const double recovery_rate_;
  const double lower_;
  const double upper_;
};

class pcds_edtl_functor : public dtl_functor {
 public:
  template <typename _ForwardIt1, typename _ForwardIt2>
  pcds_edtl_functor(_ForwardIt1 t_first, _ForwardIt1 t_last, _ForwardIt2 df_first,
                    const double recovery_rate, const double coupon, const double upfront)
      : dtl_functor(
            eddl_functor{df_first, std::next(df_first, std::distance(t_first, t_last))}.coeffs(),
            pcds_edpl_functor{t_first, t_last, df_first, recovery_rate, coupon, upfront}.coeffs()) {
  }

  template <typename _ForwardIt>
  double operator()(_ForwardIt l_first, _ForwardIt l_last) const {
    const auto l_map = [](const auto val) { return val; };

    return dl_functor::operator()(l_first, l_last, l_map);
  }
};

class cdo_edtl_functor : public dtl_functor {
 public:
  template <typename _ForwardIt1, typename _ForwardIt2>
  cdo_edtl_functor(_ForwardIt1 t_first, _ForwardIt1 t_last, _ForwardIt2 df_first,
                   const double recovery_rate, const double lower, const double upper,
                   const double coupon, const double upfront)
      : dtl_functor(
            eddl_functor{df_first, std::next(df_first, std::distance(t_first, t_last))}.coeffs(),
            cdo_edpl_functor{t_first, t_last, df_first, recovery_rate, lower, upper, coupon,
                             upfront}
                .coeffs()) {}

  template <typename _ForwardIt>
  double operator()(_ForwardIt l_first, _ForwardIt l_last) const {
    const auto l_map = [](const auto val) { return val; };

    return dl_functor::operator()(l_first, l_last, l_map);
  }
};

}  // namespace cvalr
