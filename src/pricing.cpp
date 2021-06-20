#include <Rcpp.h>
#include <cvalr.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
double Rcpp__pcds_ddl(const NumericVector &l, const NumericVector &df, const double recovery_rate) {
  using dl_functor = cvalr::pcds_ddl_functor;
  using deriv_type = typename dl_functor::deriv_type;

  return dl_functor{deriv_type{recovery_rate}, df.cbegin(), df.cend()}(l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__cdo_ddl(const NumericVector &l, const NumericVector &df, const double recovery_rate,
                     const double lower, const double upper) {
  using dl_functor = cvalr::cdo_ddl_functor;
  using deriv_type = typename dl_functor::deriv_type;

  return dl_functor{deriv_type{recovery_rate, lower, upper}, df.cbegin(), df.cend()}(l.cbegin(),
                                                                                     l.cend());
}

// [[Rcpp::export]]
double Rcpp__eddl(const NumericVector &l, const NumericVector &df) {
  using dl_functor =
      cvalr::pcds_eddl_functor;  // choice arbitrary, could also be cdo_eddl_functor
  using deriv_type = typename dl_functor::deriv_type;

  return dl_functor{deriv_type{0.}, df.cbegin(), df.cend()}(l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__pcds_dpl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                      const double recovery_rate, const double coupon, const double upfront) {
  using dl_functor = cvalr::pcds_dpl_functor;
  using deriv_type = typename dl_functor::deriv_type;

  return dl_functor{deriv_type{recovery_rate, coupon, upfront}, t.cbegin(), t.cend(), df.cbegin()}(
      l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__cdo_dpl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                     const double recovery_rate, const double lower, const double upper,
                     const double coupon, const double upfront) {
  using dl_functor = cvalr::cdo_dpl_functor;
  using deriv_type = typename dl_functor::deriv_type;
  return dl_functor{deriv_type{recovery_rate, lower, upper, coupon, upfront}, t.cbegin(), t.cend(),
                    df.cbegin()}(l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__pcds_edpl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                       const double recovery_rate, const double coupon, const double upfront) {
  using dl_functor = cvalr::pcds_edpl_functor;
  using deriv_type = typename dl_functor::deriv_type;

  return dl_functor{deriv_type(recovery_rate, coupon, upfront), t.cbegin(), t.cend(), df.cbegin()}(
      l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__cdo_edpl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                      const double recovery_rate, const double lower, const double upper,
                      const double coupon, const double upfront) {
  using dl_functor = cvalr::cdo_edpl_functor;
  using deriv_type = typename dl_functor::deriv_type;
  return dl_functor{deriv_type{recovery_rate, lower, upper, coupon, upfront}, t.cbegin(), t.cend(),
                    df.cbegin()}(l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__pcds_dtl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                      const double recovery_rate, const double coupon, const double upfront) {
  using dl_functor = cvalr::pcds_dtl_functor;
  using deriv_type = typename dl_functor::deriv_type;
  return dl_functor{deriv_type{recovery_rate, coupon, upfront}, t.cbegin(), t.cend(), df.cbegin()}(
      l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__cdo_dtl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                     const double recovery_rate, const double lower, const double upper,
                     const double coupon, const double upfront) {
  using dl_functor = cvalr::cdo_dtl_functor;
  using deriv_type = typename dl_functor::deriv_type;
  return dl_functor{deriv_type{recovery_rate, lower, upper, coupon, upfront}, t.cbegin(), t.cend(),
                    df.cbegin()}(l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__pcds_edtl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                       const double recovery_rate, const double coupon, const double upfront) {
  using dl_functor = cvalr::pcds_edtl_functor;
  using deriv_type = typename dl_functor::deriv_type;
  return dl_functor{deriv_type{recovery_rate, coupon, upfront}, t.cbegin(), t.cend(), df.cbegin()}(
      l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__cdo_edtl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                      const double recovery_rate, const double lower, const double upper,
                      const double coupon, const double upfront) {
  using dl_functor = cvalr::cdo_edtl_functor;
  using deriv_type = typename dl_functor::deriv_type;
  return dl_functor{deriv_type{recovery_rate, lower, upper, coupon, upfront}, t.cbegin(), t.cend(),
                    df.cbegin()}(l.cbegin(), l.cend());
}
