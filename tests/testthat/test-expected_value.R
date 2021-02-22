times <- seq(0, 5, by = 0.25)

dim = 50L
lambda <- 8e-2
rho <- 4e-1

seed <- 1623
n_sim <- 1e2

recovery_rate <- 0.4
lower <- 0.1
upper <- 0.2

lower_vec <- c(0, 0.1, 0.2, 0.35)
upper_vec <- c(0.1, 0.2, 0.35, 1)

test_that("`expected_value` works as expected for `CalibrationParam`", {
  parm <- AlphaStableExtMO2FParam(dim = dim, lambda = lambda, rho = rho)

  g <- function(k, recovery_rate, lower, upper) {
    pmin(pmax((1 - recovery_rate) * k - lower, 0), upper - lower)
  }
  x <- expected_value(parm, times, g,
    lower = lower, upper = upper, recovery_rate = recovery_rate,
    pd_args = list(method = "ExMarkovParam", seed = seed, n_sim = n_sim))
  expect_numeric(x, any.missing = FALSE, lower = 0, upper = 1,
    len = length(times), sorted = TRUE)
})
