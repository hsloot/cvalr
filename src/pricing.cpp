#include <Rcpp.h>
#include <cvalr.hpp>
#include <fastcvalr.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
double Rcpp__pcds_ddl(const NumericVector &l, const NumericVector &df, const double recovery_rate) {
  return cvalr::pcds_ddl_functor{df.cbegin(), df.cend(), recovery_rate}(l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__cdo_ddl(const NumericVector &l, const NumericVector &df, const double recovery_rate,
                     const double lower, const double upper) {
  return cvalr::cdo_ddl_functor{df.cbegin(), df.cend(), recovery_rate, lower, upper}(l.cbegin(),
                                                                                     l.cend());
}

// [[Rcpp::export]]
double Rcpp__eddl(const NumericVector &l, const NumericVector &df) {
  return cvalr::eddl_functor{df.cbegin(), df.cend()}(l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__pcds_dpl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                      const double recovery_rate, const double coupon, const double upfront) {
  return cvalr::pcds_dpl_functor{t.cbegin(), t.cend(), df.cbegin(), recovery_rate, coupon, upfront}(
      l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__cdo_dpl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                     const double recovery_rate, const double lower, const double upper,
                     const double coupon, const double upfront) {
  return cvalr::cdo_dpl_functor{t.cbegin(), t.cend(), df.cbegin(), recovery_rate,
                                lower,      upper,    coupon,      upfront}(l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__pcds_edpl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                       const double recovery_rate, const double coupon, const double upfront) {
  return cvalr::pcds_edpl_functor{t.cbegin(),    t.cend(), df.cbegin(),
                                  recovery_rate, coupon,   upfront}(l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__cdo_edpl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                      const double recovery_rate, const double lower, const double upper,
                      const double coupon, const double upfront) {
  return cvalr::cdo_edpl_functor{t.cbegin(), t.cend(), df.cbegin(), recovery_rate,
                                 lower,      upper,    coupon,      upfront}(l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__pcds_dtl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                      const double recovery_rate, const double coupon, const double upfront) {
  return cvalr::pcds_dtl_functor{t.cbegin(), t.cend(), df.cbegin(), recovery_rate, coupon, upfront}(
      l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__cdo_dtl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                     const double recovery_rate, const double lower, const double upper,
                     const double coupon, const double upfront) {
  return cvalr::cdo_dtl_functor{t.cbegin(), t.cend(), df.cbegin(), recovery_rate,
                                lower,      upper,    coupon,      upfront}(l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__pcds_edtl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                       const double recovery_rate, const double coupon, const double upfront) {
  return cvalr::pcds_edtl_functor{t.cbegin(),    t.cend(), df.cbegin(),
                                  recovery_rate, coupon,   upfront}(l.cbegin(), l.cend());
}

// [[Rcpp::export]]
double Rcpp__cdo_edtl(const NumericVector &l, const NumericVector &t, const NumericVector &df,
                      const double recovery_rate, const double lower, const double upper,
                      const double coupon, const double upfront) {
  return cvalr::cdo_edtl_functor{t.cbegin(), t.cend(), df.cbegin(), recovery_rate,
                                 lower,      upper,    coupon,      upfront}(l.cbegin(), l.cend());
}
