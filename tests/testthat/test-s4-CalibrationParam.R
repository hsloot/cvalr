d <- 5L
lambda <- 8e-2
alpha <- 4e-1

rho <- 3 * alpha / (4 - alpha)
tau <- alpha / (2 - alpha)


n <- 5e1L
times <- seq(0, 5, by = 0.25)

test_that("`probability_distribution` works as expected for `CalibrationParam`", {
  # HELPER START
  pfn <- function(parm, times, n, seed = NULL) {
    d <- getDimension(parm)
    if (!is.null(seed)) {
      set.seed(seed)
    }
    simulate_adcp(parm, times, n_sim = n) %>%
      purrr::array_branch(2) %>%
      purrr::cross2((0:d) / d, .) %>%
      purrr::map_dbl(~mean(.[[2]] == .[[1]])) %>%
      matrix(nrow = d+1)
  }
  # HELPER END

  parm <- AlphaStableExtMO2FParam(d, lambda, rho = rho)

  x <- probability_distribution(
    parm, times, seed = 1623, method = "CalibrationParam", sim_args = list(n_sim = n))
  expect_matrix(
    x, mode = "numeric", any.missing = FALSE, nrows = d+1L, ncols = length(times))
  expect_numeric(
    x, lower = 0, upper = 1, finite = TRUE)
  purrr::array_branch(x, 2) %>%
    purrr::map_dbl(sum) %>%
    expect_equal(rep(1, length(times)))

  y <- pfn(parm, times, n, seed = 1623)
  expect_equal(x, y)

  y <- probability_distribution(
    parm, times, seed = 1623, method = "CalibrationParam",
    sim_args = list(n_sim = n, method = "ExMOParam"))
  expect_equal(x, y)

  x <- probability_distribution(
    parm, times, seed = 1623, method = "CalibrationParam",
    sim_args = list(n_sim = n, method = "ExMarkovParam"))
  y <- pfn(as(parm, "ExMarkovParam"), times, n, seed = 1623)
  expect_equal(x, y)
})


g_pcds <- function(x, recovery_rate) {
  (1 - recovery_rate) * x
}
g_cdo <- function(x, recovery_rate, lower, upper) {
  pmin(pmax((1 - recovery_rate) * x - lower, 0), upper - lower)
}

recovery_rate <- 4e-1
lower <- 0.15
upper <- 0.25

# TODO: Use a different `g` here (maybe nth-to-default).
test_that("`expected_value` works as expected for `CalibrationParam`", {
  # HELPER START
  efn1 <- function(t, g, parm) {
    mu <- sapply((0L:d) / d, g, recovery_rate = recovery_rate, lower = lower, upper = upper)
    probs <- probability_distribution(parm, t)

    as.vector(t(probs) %*% mu)
  }
  efn2 <- function(t, g, parm, n, seed) {
    mu <- sapply((0L:d) / d, g, recovery_rate = recovery_rate, lower = lower, upper = upper)
    probs <- probability_distribution(
      parm, t, method = "CalibrationParam", seed = seed, sim_args = list(n_sim = n))

    as.vector(t(probs) %*% mu)
  }
  # HELPER END

  parm <- ExtGaussian2FParam(d, lambda, rho = rho)

  # using the default
  x <- expected_value(
    parm, times, g_cdo, lower = lower, upper = upper, recovery_rate = recovery_rate)
  expect_numeric(
    x, any.missing = FALSE, lower = 0, upper = 1, len = length(times), sorted = TRUE)

  y <- efn1(times, g_cdo, parm)
  expect_equal(x, y)

  # using MC-simulation
  x <- expected_value(
    parm, times, g_cdo, lower = lower, upper = upper, recovery_rate = recovery_rate,
    pd_args = list(method = "CalibrationParam", seed = 1623, sim_args = list(n_sim = n)))
  expect_numeric(
    x, any.missing = FALSE, lower = 0, upper = 1, len = length(times), sorted = TRUE)

  y <- efn2(times, g_cdo, parm, n, 1623)
  expect_equal(x, y)
})


test_that("`expected_pcds_loss` works as expected for `CalibrationParam`", {
  parm <- ExtGaussian2FParam(dim = d, lambda = lambda, rho = rho)

  x <- expected_pcds_loss(
    parm, times, recovery_rate = recovery_rate, method = "CalibrationParam")
  expect_numeric(x, lower = 0, upper = 1, any.missing = FALSE, len = length(times), sorted = TRUE)
  expect_equal(x, expected_value(parm, times, g_pcds, recovery_rate = recovery_rate))
  expect_equal(x, (1 - recovery_rate) * pexp(times, rate = lambda))
})


test_that("`expected_cdo_loss` works as expected for `CalibrationParam`", {
  parm <- ExtGaussian2FParam(dim = d, lambda = lambda, rho = rho)

  x <- expected_cdo_loss(
    parm, times, recovery_rate = recovery_rate, lower = lower, upper = upper,
    method = "CalibrationParam")
  expect_numeric(x, lower = 0, upper = 1, any.missing = FALSE, len = length(times), sorted = TRUE)
  expect_equal(
    x, expected_value(
      parm, times, g_cdo, recovery_rate = recovery_rate, lower = lower, upper = upper))
})


discount_factors <- exp(-times * 50e-4)
coupon <- 8e-2
upfront <- 1e-2

test_that("`expected_pcds_equation` works as expected for `CalibrationParam`", {
  parm <- ExtGaussian2FParam(dim = d, lambda = lambda, rho = rho)

  x <- expected_pcds_equation(
    parm, times, discount_factors, recovery_rate, coupon, upfront)
  expect_numeric(x, finite = TRUE, any.missing = FALSE, len = 1)
  y <- portfolio_cds_equation(
    expected_losses = expected_pcds_loss(parm, times, recovery_rate),
    times = times, discount_factors = discount_factors, recovery_rate = recovery_rate,
    coupon = coupon, upfront = upfront)
  expect_equal(x, y)
})


lower_vec <- c(0, 0.1, 0.2, 0.35)
upper_vec <- c(0.1, 0.2, 0.35, 1)
coupon_vec <- c(0.05, 0.05, 0.05, 0.015)
upfront_vec <- c(0.7, 0.4, 0.1, 0)

test_that("`expected_cdo_equation` works as expected for `CalibrationParam`", {
  parm <- ExtGaussian2FParam(dim = d, lambda = lambda, rho = rho)

  x <- expected_cdo_equation(
    parm, times, discount_factors, recovery_rate, lower_vec, upper_vec, coupon_vec, upfront_vec)
  expect_numeric(x, finite = TRUE, any.missing = FALSE, len = length(lower_vec))

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
