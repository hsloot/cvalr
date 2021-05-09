times <- seq(0, 5, by = 0.25)

dim <- 50L
lambda <- 8e-2
rho <- 4e-1

seed <- 1623
n_sim <- 1e2

recovery_rate <- 0.4
lower <- 0.1
upper <- 0.2

lower_vec <- c(0, 0.1, 0.2, 0.35)
upper_vec <- c(0.1, 0.2, 0.35, 1)

spread <- 8e-2
coupon_vec <- c(0.05, 0.05, 0.05, 0.015)
upfront_vec <- c(0.7, 0.4, 0.1, 0)

discount_factors <- rep(1, length(times))

g_pcds <- function(k, recovery_rate) {
  (1 - recovery_rate) * k
}

g_cdo <- function(k, recovery_rate, lower, upper) {
  pmin(pmax((1 - recovery_rate) * k - lower, 0), upper - lower)
}

test_that("`expected_value` works as expected for `CalibrationParam`", {
  parm <- AlphaStableExtMO2FParam(dim = dim, lambda = lambda, rho = rho)

  x <- expected_value(parm, times, g_cdo,
    lower = lower, upper = upper, recovery_rate = recovery_rate,
    pd_args = list(method = "CalibrationParam", seed = seed, n_sim = n_sim))
  expect_numeric(x, any.missing = FALSE, lower = 0, upper = 1,
    len = length(times), sorted = TRUE)

  mu <- sapply((0:dim) / dim, g_cdo,
    recovery_rate = recovery_rate, lower = lower, upper = upper)
  probs <- probability_distribution(parm, times, method = "CalibrationParam",
    seed = seed, n_sim = n_sim)

  expect_equal(x, as.vector(t(probs) %*% mu))
})

test_that("`expected_pcds_loss` works as intended for `CalibrationParam`", {
  parm <- ExtGaussian2FParam(dim = dim, lambda = lambda, rho = rho)

  x <- expected_pcds_loss(parm, times, recovery_rate = recovery_rate,
    method = "CalibrationParam")
  expect_numeric(x, any.missing = FALSE, lower = 0, upper = 1,
    len = length(times), sorted = TRUE)

  expect_equal(x, expected_value(parm, times, g_pcds, recovery_rate = recovery_rate))

  mu <- sapply((0:dim) / dim, g_pcds, recovery_rate = recovery_rate)
  probs <- probability_distribution(parm, times)
  expect_equal(x, as.vector(t(probs) %*% mu))
})

test_that("`expected_pcds_loss` works as intended for `ExtMO2FParam", {
  parm <- AlphaStableExtMO2FParam(dim = dim, lambda = lambda, rho = rho)

  x <- expected_pcds_loss(parm, times, recovery_rate = recovery_rate)
  expect_numeric(x, any.missing = FALSE, lower = 0, upper = 1,
    len = length(times), sorted = TRUE)
  expect_equal(x, (1 - recovery_rate) * pexp(times, rate = parm@lambda))
})

test_that("`expected_pcds_loss` works as intended for `ExtArch2FParam", {
  parm <- FrankExtArch2FParam(dim = dim, lambda = lambda, rho = rho)

  x <- expected_pcds_loss(parm, times, recovery_rate = recovery_rate)
  expect_numeric(x, any.missing = FALSE, lower = 0, upper = 1,
    len = length(times), sorted = TRUE)
  expect_equal(x, (1 - recovery_rate) * pexp(times, rate = parm@lambda))
})

test_that("`expected_cdo_loss` works as intended for `CalibrationParam`", {
  parm <- ExtGaussian2FParam(dim = dim, lambda = lambda, rho = rho)

  x <- expected_cdo_loss(parm, times, recovery_rate = recovery_rate,
    lower = lower, upper = upper, method = "CalibrationParam")
  expect_numeric(x, any.missing = FALSE, lower = 0, upper = 1,
    len = length(times), sorted = TRUE)

  expect_equal(x, expected_value(parm, times, g_cdo,
    recovery_rate = recovery_rate, lower = lower, upper = upper))

  mu <- sapply((0:dim) / dim, g_cdo,
    recovery_rate = recovery_rate, lower = lower, upper = upper)
  probs <- probability_distribution(parm, times)
  expect_equal(x, as.vector(t(probs) %*% mu))
})

test_that("`expected_pcds_equation` works as expected for `CalibrationParam`", {
  parm <- ExtGaussian2FParam(dim = dim, lambda = lambda, rho = rho)

  x <- expected_pcds_equation(parm, times, discount_factors,
    recovery_rate, spread, 0)
  expect_numeric(x, any.missing = FALSE, len = 1)

  expected_losses <- expected_pcds_loss(parm, times, recovery_rate)
  expect_equal(x,
    portfolio_cds_equation(
      expected_losses = expected_losses, times = times,
      discount_factors = discount_factors, recovery_rate = recovery_rate,
      coupon = spread, upfront = 0))
})

test_that("`expected_cdo_equation` works as expected for `CalibrationParam`", {
  parm <- ExtGaussian2FParam(dim = dim, lambda = lambda, rho = rho)

  x <- expected_cdo_equation(parm, times, discount_factors,
    recovery_rate, lower_vec, upper_vec, coupon_vec, upfront_vec)
  expect_numeric(x, any.missing = FALSE, len = 4)

  for (i in seq_along(x)) {
    expected_losses <- expected_cdo_loss(parm, times, recovery_rate,
      lower_vec[i], upper_vec[i])
    expect_equal(x[i],
      cdo_equation(expected_losses = expected_losses, times = times,
        discount_factors = discount_factors,
        lower = lower_vec[i], upper = upper_vec[i],
        coupon = coupon_vec[i], upfront = upfront_vec[i]))
  }
})
