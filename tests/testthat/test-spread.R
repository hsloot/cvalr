lambda <- 0.1
rate <- 0.005
recovery_rate <- 0.4
times <- seq(0, 2, by = 0.25)
expected_default_count <- 1 - exp(-lambda * times)
discount_factors <- exp(-rate * times)

test_that("eddl is calculated correctly", {
  expected_losses <- (1 - recovery_rate) * expected_default_count
  expected_nominals <- 1 - expected_default_count
  expect_equal(
    portfolio_cds_spread(
      expected_default_count, times, discount_factors, recovery_rate
    ),
    sum(discount_factors[-1] * diff(expected_losses)) /
      sum(
        discount_factors[-1] * diff(times) *
          0.5 * (expected_nominals[-1] + expected_nominals[-9])
      )
  )
})

test_that("upfront payment is calculated correctly", {
  spread <- 0.05
  lower <- 0.03
  upper <- 0.1
  expected_losses <- pmin(
    pmax(
      (1 - recovery_rate) * expected_default_count - lower,
      0
    ),
    upper - lower
  )
  expected_nominals <- upper - lower - expected_losses
  expect_equal(
    upfront_payment(
      expected_default_count, times, discount_factors, lower, upper, spread, recovery_rate
    ),
    (sum(discount_factors[-1] * diff(expected_losses)) - spread *
      sum(
        discount_factors[-1] * diff(times) *
          0.5 * (expected_nominals[-1] + expected_nominals[-9])
      )) / (upper - lower)
  )
})

test_that("upfront spread is calculated correctly", {
  lower <- 0.03
  upper <- 0.1
  expected_losses <- pmin(
    pmax(
      (1 - recovery_rate) * expected_default_count - lower,
      0
    ),
    upper - lower
  )
  expected_nominals <- upper - lower - expected_losses
  expect_equal(
    upfront_spread(
      expected_default_count, times, discount_factors, lower, upper, recovery_rate
    ),
    sum(discount_factors[-1] * diff(expected_losses)) /
      sum(
        discount_factors[-1] * diff(times) *
          0.5 * (expected_nominals[-1] + expected_nominals[-9])
      )
  )
})
