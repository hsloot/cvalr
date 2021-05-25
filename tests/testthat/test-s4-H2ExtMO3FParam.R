partition <- list(1L:2L, 3L:6L, 7L:8L)
composition <- purrr::map_int(partition, length)
d <- sum(composition)
lambda <- 8e-2
alpha <- c(3e-1, 6e-1)

test_that("`expected_pcds_loss` works as expected for `H2ExtMO3FParam", {
  # HELPER START
  epcdslfn <- function(parm, times, recovery_rate) {
    qassert(times, "N+[0,)")
    qassert(recovery_rate, "N1[0,1]")

    (1 - recovery_rate) * pexp(times, rate = getLambda(parm))
  }
  # HELPER END

  parm <- AlphaStableH2ExtMO3FParam(composition = composition, lambda = lambda, alpha = alpha)
  times <- seq(25e-2, 5L, by = 25e-2)
  recovery_rate <- 0.4

  x <- expected_pcds_loss(parm, times, recovery_rate = recovery_rate)
  expect_numeric(x, any.missing = FALSE, lower = 0, upper = 1,
    len = length(times), sorted = TRUE)
  expect_equal(x, epcdslfn(parm, times, recovery_rate))
})
