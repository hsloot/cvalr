d <- 75L
lambda <- 0.1

recovery_rate <- 0.4
lower <- 0.1
upper <- 0.2
coupon <- 20e-2
upfront <- -5e-2


times <- seq(25e-2, 5, by = 0.25)
df <- exp(-(-50e-4) * times)
l <- drop(Rcpp__dt2adcp(matrix(rexp(d, lambda), nrow = 1L), times))

test_that("Portfolio CDS DDL is calculated correctly", {
  expect_equal(Rcpp__pcds_ddl(l, df, recovery_rate),
               test__pcds_ddl(l, df, recovery_rate))
})

test_that("CDO DDL is calculated correctly", {
  expect_equal(Rcpp__cdo_ddl(l, df, recovery_rate, lower, upper),
               test__cdo_ddl(l, df, recovery_rate, lower, upper))
})

test_that("EDDL is calculated correctly", {
  el <- pmin(pmax((1- recovery_rate) * l - lower, 0), upper - lower)

  expect_equal(Rcpp__eddl(el, df),
               test__eddl(el, df))
})

test_that("Portfolio CDS DPL is calculated correctly", {
  expect_equal(Rcpp__pcds_dpl(l, times, df, recovery_rate, coupon, upfront),
               test__pcds_dpl(l, times, df, recovery_rate, coupon, upfront))
})

test_that("CDO DPL is calculated correctly", {
  expect_equal(Rcpp__cdo_dpl(l, times, df, recovery_rate, lower, upper, coupon, upfront),
               test__cdo_dpl(l, times, df, recovery_rate, lower, upper, coupon, upfront))
})

test_that("Portfolio CDS EDPL is calculated correctly", {
  expect_equal(Rcpp__pcds_edpl(l, times, df, recovery_rate, coupon, upfront),
               test__pcds_edpl(l, times, df, recovery_rate, coupon, upfront))
})

test_that("CDO EDPL is calculated correctly", {
  expect_equal(Rcpp__cdo_edpl(l, times, df, recovery_rate, lower, upper, coupon, upfront),
               test__cdo_edpl(l, times, df, recovery_rate, lower, upper, coupon, upfront))
})

test_that("Portfolio CDS DTL is calculated correctly", {
  expect_equal(Rcpp__pcds_dtl(l, times, df, recovery_rate, coupon, upfront),
               test__pcds_dtl(l, times, df, recovery_rate, coupon, upfront))
})

test_that("CDO DTL is calculated correctly", {
  expect_equal(Rcpp__cdo_dtl(l, times, df, recovery_rate, lower, upper, coupon, upfront),
               test__cdo_dtl(l, times, df, recovery_rate, lower, upper, coupon, upfront))
})

test_that("Portfolio CDS EDTL is calculated correctly", {
  expect_equal(Rcpp__pcds_edtl(l, times, df, recovery_rate, coupon, upfront),
               test__pcds_edtl(l, times, df, recovery_rate, coupon, upfront))
})

test_that("CDO EDTL is calculated correctly", {
  expect_equal(Rcpp__cdo_edtl(l, times, df, recovery_rate, lower, upper, coupon, upfront),
               test__cdo_edtl(l, times, df, recovery_rate, lower, upper, coupon, upfront))
})
