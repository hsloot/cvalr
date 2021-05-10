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



test_that("`expected_pcds_loss` works as intended for `ExtArch2FParam", {
  parm <- FrankExtArch2FParam(dim = dim, lambda = lambda, rho = rho)

  x <- expected_pcds_loss(parm, times, recovery_rate = recovery_rate)
  expect_numeric(x, any.missing = FALSE, lower = 0, upper = 1,
    len = length(times), sorted = TRUE)
  expect_equal(x, (1 - recovery_rate) * pexp(times, rate = parm@lambda))
})
