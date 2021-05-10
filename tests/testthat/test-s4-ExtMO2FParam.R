d <- 5L
lambda <- 8e-2
rho <- 4e-1


times <- seq(0, 5, by = 0.25)
recovery_rate <- 4e-1

test_that("`expected_pcds_loss` works as expected for `ExtMO2FParam", {
  parm <- AlphaStableExtMO2FParam(dim = d, lambda = lambda, rho = rho)

  x <- expected_pcds_loss(parm, times, recovery_rate = recovery_rate)
  expect_numeric(
    x, any.missing = FALSE, lower = 0, upper = 1, len = length(times), sorted = TRUE)
  expect_equal(x, (1 - recovery_rate) * pexp(times, rate = lambda))
})
