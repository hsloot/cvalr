lambda <- 0.1
rate <- 0.005
recovery_rate <- 0.4
times <- seq(0, 2, by = 0.25)
discount_factors <- exp(-rate * times)

expected_losses <- cds_expected_loss(times, lambda, recovery_rate)

test_that("Portfolio CDS coupon is calculated correctly", {
  expected_nominals <- 1 - expected_losses / (1 - recovery_rate)

  eddl <- sum(discount_factors[-1] * diff(expected_losses))
  edpl1 <- sum(
    discount_factors[-1] * diff(times) *
      0.5 * (expected_nominals[-1] + expected_nominals[-9])
  )

  expect_equal(
    portfolio_cds_coupon(
      expected_losses, times, discount_factors, recovery_rate
    ),
    eddl / edpl1
  )
})

test_that("Portfolio CDS upfront is calculated correctly", {
  coupon <- 0.05
  expected_nominals <- 1 - expected_losses / (1 - recovery_rate)

  eddl <- sum(discount_factors[-1] * diff(expected_losses))
  edpl1 <- sum(
    discount_factors[-1] * diff(times) *
      0.5 * (expected_nominals[-1] + expected_nominals[-9])
  )

  expect_equal(
    portfolio_cds_upfront(
      expected_losses, times, discount_factors, recovery_rate, coupon
    ),
    eddl - (coupon * edpl1)
  )
})

test_that("Portfolio CDS equation is calculated correctly", {
  coupon <- 0.05
  upfront <- portfolio_cds_upfront(
    expected_losses, times, discount_factors, recovery_rate, coupon
  )

  expect_equal(
    portfolio_cds_equation(
      expected_losses, times, discount_factors, recovery_rate, coupon, upfront
    ),
    0
  )
})

lower <- 0.03
upper <- 0.1
rho <- 0.1
expected_losses <- cdo_tranche_expected_loss_gaussian(times, lambda, rho, lower, upper, recovery_rate)
expected_nominals <- upper - lower - expected_losses

test_that("CDO upfront payment is calculated correctly", {
  coupon <- 0.05

  eddl <- sum(discount_factors[-1] * diff(expected_losses))
  edpl1 <- sum(
    discount_factors[-1] * diff(times) *
      0.5 * (expected_nominals[-1] + expected_nominals[-9])
  )

  expect_equal(
    cdo_upfront(
      expected_losses, times, discount_factors, lower, upper, coupon
    ),
    (eddl - coupon * edpl1) / (upper - lower))
})

test_that("CDO coupon is calculated correctly", {

  eddl <- sum(discount_factors[-1] * diff(expected_losses))
  edpl1 <- sum(
    discount_factors[-1] * diff(times) *
      0.5 * (expected_nominals[-1] + expected_nominals[-9])
  )

  expect_equal(
    cdo_coupon(
      expected_losses, times, discount_factors, lower, upper
    ),
    eddl / edpl1
  )
})

test_that("CDO equation is calculated correctly", {
  coupon <- 0.05
  upfront <- cdo_upfront(
    expected_losses, times, discount_factors, lower, upper, coupon
  )
  expect_equal(
    cdo_equation(
      expected_losses, times, discount_factors, lower, upper, coupon, upfront
    ),
    0
  )
})
