test_that("eddl is calculated correctly", {
  lambda <- 0.1
  rate <- 0.005
  recovery_rate <- 0.4
  times <- seq(0, 2, by = 0.25)
  expected_default_count <- 1 - exp(-lambda * times)
  expected_nominals <- 1 - expected_default_count
  expected_losses <- (1 - recovery_rate) * expected_default_count
  discount_factors <- exp(-rate * times)
  expect_equal(
    portfolio_cds_spread(
      expected_default_count, times, discount_factors, recovery_rate),
    sum(discount_factors[-1] * diff(expected_losses)) /
      sum(
        discount_factors[-1] * diff(times) *
          0.5 * (expected_nominals[-1] + expected_nominals[-9]))
  )
})
