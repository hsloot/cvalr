set.seed(1623)

d <- 5L
lambda <- 8e-2
alpha <- 4e-1

rho <- 3 * alpha / (4 - alpha)
tau <- alpha / (2 - alpha)

n_sim <- 5e1L
times <- seq(25e-2, 5, by = 25e-2)
discount_factors <- exp(-times * 50e-4)

recovery_rate <- 4e-1

test_that("`expected_pcds_equation` works as expected for `CalibrationParam`", {
  coupon <- 1e-1
  upfront <- -1e-2
  parm <- as(ArmageddonExtMO2FParam(d, lambda, rho = rho), "ExMOParam")

  # using default
  x <- expected_pcds_equation(
    parm, times, discount_factors, recovery_rate, coupon, upfront,
    method = "prob")
  expect_numeric(x, finite = TRUE, any.missing = FALSE, len = 1L)
  y <- test__expected_pcds_equation__prob(
    parm, times, discount_factors, recovery_rate, coupon, upfront)
  expect_equal(x, y)

  # using MC-simulation
  set.seed(1623)
  x <- expected_pcds_equation(
    parm, times, discount_factors, recovery_rate, coupon, upfront,
    method = "mc", n_sim = n_sim)
  expect_numeric(x, finite = TRUE, any.missing = FALSE, len = 1L)
  set.seed(1623)
  y <- test__expected_pcds_equation__mc(
    parm, times, discount_factors, recovery_rate, coupon, upfront,
    n_sim = n_sim)
  expect_equal(x, y)
})

test_that("`expected_cdo_equation` works as expected for `CalibrationParam`", {
  lower <- c(0, 0.1, 0.2, 0.35)
  upper <- c(0.1, 0.2, 0.35, 1)
  coupon <- c(rep(5e-2, 3), 0)
  upfront <- c(8e-1, 5e-1, 1e-1, 0)
  parm <- as(ArmageddonExtMO2FParam(d, lambda, rho = rho), "ExMOParam")

  # using default
  x <- expected_cdo_equation(
    parm, times, discount_factors, recovery_rate, lower, upper, coupon, upfront,
    method = "prob")
  expect_numeric(x, finite = TRUE, any.missing = FALSE, len = 4L)
  y <- test__expected_cdo_equation__prob(
    parm, times, discount_factors, recovery_rate, lower, upper, coupon, upfront)
  expect_equal(x, y)

  # using MC-simulation
  set.seed(1623)
  x <- expected_cdo_equation(
    parm, times, discount_factors, recovery_rate, lower, upper, coupon, upfront,
    method = "mc", n_sim = n_sim)
  expect_numeric(x, finite = TRUE, any.missing = FALSE, len = 4L)
  set.seed(1623)
  y <- test__expected_cdo_equation__mc(
    parm, times, discount_factors, recovery_rate, lower, upper, coupon, upfront,
    n_sim = n_sim)
  expect_equal(x, y)
})
