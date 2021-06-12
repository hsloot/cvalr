partition <- list(1L:2L, 3L:6L, 7L:8L)
composition <- purrr::map_int(partition, length)
d <- sum(composition)
lambda <- 8e-2
tau <- c(3e-1, 6e-1)

test_that("`expected_pcds_equation` works as expected for `H2ExtArch3FParam", {
  times <- seq(25e-2, 5, by = 25e-2)
  discount_factors <- rep(1, length(times))
  recovery_rate <- 4e-1
  coupon <- 1e-1
  upfront <- -1e-2
  parm <- H2ExtArch3FParam(composition = composition, lambda = lambda, tau = tau)

  # using default
  x <- expected_pcds_equation(parm, times, discount_factors, recovery_rate, coupon, upfront)
  expect_numeric(x, finite = TRUE, any.missing = FALSE, len = 1L)
  y <- test__expected_pcds_equation__default(
    parm, times, discount_factors, recovery_rate, coupon, upfront)
  expect_equal(x, y)
})
