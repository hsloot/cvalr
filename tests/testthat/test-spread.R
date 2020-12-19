test_that("eddl is calculated correctly", {
  lambda <- 0.1
  rate <- 0.005
  recovery_rate <- 0.4
  times <- seq(0, 2, by = 0.25)
  expected_nominals <- exp(-lambda * times)
  expected_losses <- (1 - recovery_rate) * (1 - expected_nominals)
  discount_factors <- exp(-rate * times)
  expect_equal(
    portfolio_cds_spread(
      expected_losses, expected_nominals, times, discount_factors),
    sum(discount_factors[-1] * diff(expected_losses)) /
      sum(
        discount_factors[-1] * diff(times) *
          0.5 * (expected_nominals[-1] + expected_nominals[-9]))
  )
})
