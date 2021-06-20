#pragma once

#include "cvalr/conversion.hpp"
#include "cvalr/derivative.hpp"
#include "cvalr/payments.hpp"

namespace cvalr {

  using pcds_derivative = derivative::pcds_derivative;
  using cdo_derivative = derivative::cdo_derivative;
  using expected_pcds_derivative = derivative::expectation_facade<derivative::pcds_derivative>;
  using expected_cdo_derivative = derivative::expectation_facade<derivative::cdo_derivative>;

  using pcds_ddl_functor = payments::ddl_functor<pcds_derivative>;
  using cdo_ddl_functor = payments::ddl_functor<cdo_derivative>;

  using pcds_eddl_functor = payments::ddl_functor<expected_pcds_derivative>;
  using cdo_eddl_functor = payments::ddl_functor<expected_cdo_derivative>;

  using pcds_dpl_functor = payments::dpl_functor<pcds_derivative>;
  using cdo_dpl_functor = payments::dpl_functor<cdo_derivative>;

  using pcds_edpl_functor = payments::dpl_functor<expected_pcds_derivative>;
  using cdo_edpl_functor = payments::dpl_functor<expected_cdo_derivative>;

  using pcds_dtl_functor = payments::dtl_functor<pcds_derivative>;
  using cdo_dtl_functor = payments::dtl_functor<cdo_derivative>;

  using pcds_edtl_functor = payments::dtl_functor<expected_pcds_derivative>;
  using cdo_edtl_functor = payments::dtl_functor<expected_cdo_derivative>;

  using conversion::dt2adcp;

}  // namespace cvalr
