## #### Global settings ####
n_sim <- 1e3

n_tests <- 2L * 8L + 9L

level <- 5e-2
level <- level / n_tests

set.seed(1632)


# see https://en.wikipedia.org/wiki/Fisher_transformation
# and https://stats.stackexchange.com/a/506367
expect_cor_not_rejected <- function(x, corr, level) {
  n <- if (is.vector(x)) length(x) else nrow(x)
  theta <- copula::P2p(corr)
  eta <- atanh(theta)
  thetahat <- copula::P2p(cor(x = x, method = "spearman"))
  etahat <- atanh(thetahat)
  c <- sqrt(1 + theta^2/2)
  b <- 3
  z <- pnorm(etahat, mean = eta, sd = c / sqrt(n - b))

  testthat::expect(
    ok = max(abs(z)) <= 1 - level / (2 * length(z)),
    failure_message = sprintf(
      "Null hypothesis rejected %e < %e",
      max(abs(z)), 1 - level / (2 * length(z))))
}

expect_dist_not_rejected <- function(x, lambda, level) {
  p <- suppressWarnings(apply(x, 2L, function(x) ks.test(x, "pexp", rate = lambda)$p.value))
  testthat::expect(
    ok = max(p) >= level / length(p),
    failure_message = sprintf(
      "Null hypothesis rejected %e < %e", max(p), level / length(p)
    ))
}

composition <- c(28, 24, 16, 4, 3)
d <- sum(composition)
lambda <- 8e-2
rho <- 4e-1
tau <- 4e-1
corr <- copula::p2P(rho, d = d)

test_that("`simulate_dt` passes marginal hypothesis test for ext. params", {
  expect_dist_not_rejected(
    simulate_dt(CuadrasAugeExtMO2FParam(dim = d, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)
  expect_dist_not_rejected(
    simulate_dt(AlphaStableExtMO2FParam(dim = d, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)
  expect_dist_not_rejected(
    simulate_dt(ExponentialExtMO2FParam(dim = d, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)
  expect_dist_not_rejected(
    simulate_dt(PoissonExtMO2FParam(dim = d, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)

  expect_dist_not_rejected(
    simulate_dt(ExtGaussian2FParam(dim = d, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)

  expect_dist_not_rejected(
    simulate_dt(ClaytonExtArch2FParam(dim = d, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)
  expect_dist_not_rejected(
    simulate_dt(FrankExtArch2FParam(dim = d, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)
  expect_dist_not_rejected(
    simulate_dt(GumbelExtArch2FParam(dim = d, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)
  expect_dist_not_rejected(
    simulate_dt(JoeExtArch2FParam(dim = d, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)
})

test_that("`simulate_dt` spearman correlation is in appr. CI for Ext. Param", {
  expect_cor_not_rejected(
    simulate_dt(CuadrasAugeExtMO2FParam(dim = d, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
  expect_cor_not_rejected(
    simulate_dt(AlphaStableExtMO2FParam(dim = d, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
  expect_cor_not_rejected(
    simulate_dt(PoissonExtMO2FParam(dim = d, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
  expect_cor_not_rejected(
    simulate_dt(ExponentialExtMO2FParam(dim = d, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)

  expect_cor_not_rejected(
    simulate_dt(ExtGaussian2FParam(dim = d, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)

  expect_cor_not_rejected(
    simulate_dt(ClaytonExtArch2FParam(dim = d, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
  expect_cor_not_rejected(
    simulate_dt(FrankExtArch2FParam(dim = d, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
  expect_cor_not_rejected(
    simulate_dt(GumbelExtArch2FParam(dim = d, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
})

fraction <- 0.4
lambda <- 10e-2
tau <- c(35e-2, 45e-2)
rho <- c(35e-2, 45e-2)
partition <- getPartition(CuadrasAugeH2ExtMO3FParam(
  composition = composition, lambda = lambda, rho = rho))

corr <- copula::p2P(rho[[1]], d = d)
walk(partition, ~{
  corr[.x, .x] <<- copula::p2P(rho[[2]], d = length(.x))
})

test_that("`simulate_dt` passes marginal hypothesis test for H2-ext. params", {
  expect_dist_not_rejected(
    simulate_dt(CuadrasAugeH2ExtMO3FParam(composition = composition, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)
  expect_dist_not_rejected(
    simulate_dt(AlphaStableH2ExtMO3FParam(composition = composition, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)
  expect_dist_not_rejected(
    simulate_dt(ExponentialH2ExtMO3FParam(composition = composition, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)
  expect_dist_not_rejected(
    simulate_dt(PoissonH2ExtMO3FParam(composition = composition, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)

  expect_dist_not_rejected(
    simulate_dt(H2ExtGaussian3FParam(composition = composition, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)

  expect_dist_not_rejected(
    simulate_dt(ClaytonH2ExtArch3FParam(composition = composition, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)
  expect_dist_not_rejected(
    simulate_dt(FrankH2ExtArch3FParam(composition = composition, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)
  expect_dist_not_rejected(
    simulate_dt(GumbelH2ExtArch3FParam(composition = composition, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)
  expect_dist_not_rejected(
    simulate_dt(JoeH2ExtArch3FParam(composition = composition, lambda = lambda, tau = tau), n_sim = n_sim),
    lambda, level = level)
})

test_that("`simulate_dt` spearman correlation is in appr. CI for H2-Ext. Param", {
  expect_cor_not_rejected(
    simulate_dt(CuadrasAugeH2ExtMO3FParam(
      composition = composition, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
  expect_cor_not_rejected(
    simulate_dt(AlphaStableH2ExtMO3FParam(
      composition = composition, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
  expect_cor_not_rejected(
    simulate_dt(PoissonH2ExtMO3FParam(
      composition = composition, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
  expect_cor_not_rejected(
    simulate_dt(ExponentialH2ExtMO3FParam(
      composition = composition, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)

  expect_cor_not_rejected(
    simulate_dt(H2ExtGaussian3FParam(
      composition = composition, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)

  expect_cor_not_rejected(
    simulate_dt(ClaytonH2ExtArch3FParam(
      composition = composition, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
  expect_cor_not_rejected(
    simulate_dt(FrankH2ExtArch3FParam(
      composition = composition, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
  expect_cor_not_rejected(
    simulate_dt(GumbelH2ExtArch3FParam(
      composition = composition, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
})
