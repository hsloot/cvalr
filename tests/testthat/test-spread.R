lambda <- 0.1
rate <- 0.005
recovery_rate <- 0.4
times <- seq(0, 2, by = 0.25)
discount_factors <- exp(-rate * times)

test_that("portfolio_cds_spread is calculated correctly", {
  expected_losses <- cds_expected_loss(times, lambda, recovery_rate)
  expected_nominals <- 1 - expected_losses / (1 - recovery_rate)
  expect_equal(
    portfolio_cds_spread(
      expected_losses, times, discount_factors, recovery_rate
    ),
    sum(discount_factors[-1] * diff(expected_losses)) /
      sum(
        discount_factors[-1] * diff(times) *
          0.5 * (expected_nominals[-1] + expected_nominals[-9])
      )
  )
})

lower <- 0.03
upper <- 0.1
rho <- 0.1
expected_losses <- cdo_tranche_expected_loss_gaussian(times, lambda, rho, lower, upper, recovery_rate)
expected_nominals <- upper - lower - expected_losses

test_that("upfront payment is calculated correctly", {
  spread <- 0.05
  expect_equal(
    upfront_payment(
      expected_losses, times, discount_factors, lower, upper, spread
    ),
    (sum(discount_factors[-1] * diff(expected_losses)) - spread *
      sum(
        discount_factors[-1] * diff(times) *
          0.5 * (expected_nominals[-1] + expected_nominals[-9])
      )) / (upper - lower)
  )
})

test_that("upfront spread is calculated correctly", {
  expect_equal(
    upfront_spread(
      expected_losses, times, discount_factors, lower, upper
    ),
    sum(discount_factors[-1] * diff(expected_losses)) /
      sum(
        discount_factors[-1] * diff(times) *
          0.5 * (expected_nominals[-1] + expected_nominals[-9])
      )
  )
})
