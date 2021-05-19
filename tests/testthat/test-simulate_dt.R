## #### Global settings ####
n_sim <- 1e3

n_tests <- 16L

level <- 5e-2
level <- level / n_tests

set.seed(1623)


# see https://en.wikipedia.org/wiki/Fisher_transformation
# and https://stats.stackexchange.com/a/506367
cor_ext <- function(x, y = NULL, use = "everything", method = "spearman", lambda = 0.95) {
  method <- match.arg(method)
  n <- if (is.vector(x)) length(x) else nrow(x)
  rs <- copula::P2p(cor(x = x, y = y, use = use, method = method))
  copula::p2P(tanh(qnorm(lambda, mean = atanh(rs), sd = sqrt((1 + rs^2) / (n-3)))))
}

cor_lower <- function(x, y = NULL, use = "everything", method = "spearman", level = 0.05) {
  cor_ext(x, y, use, method, lambda = level  / 2)
}

cor_upper <- function(x, y = NULL, use = "everything", method = "spearman", level = 0.05) {
  cor_ext(x, y, use, method, lambda = 1 - level  / 2)
}

expect_cor_not_rejected <- function(x, corr, level) {
  cp_true <- copula::P2p(corr)
  bonf_level <- level / length(cp_true)
  cp_lower <- copula::P2p(cor_lower(x, method = "spearman", level = bonf_level))
  cp_upper <- copula::P2p(cor_upper(x, method = "spearman", level = bonf_level))

  testthat::expect(
    ok = all(cp_lower <= cp_true & cp_true <= cp_upper),
    failure_message = sprintf(
      "True correlation values not in appr. Bonferroni-CI with level %.4f",
      level
    )
  )
}

d <- 8L
lambda <- 8e-2
rho <- 4e-1

corr <- copula::p2P(rho, d = d)

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
partition <- list(1:2, 3:6, 7:8)
composition <- purrr::map_int(partition, length)
d <- sum(composition)
lambda <- 8e-2
rho <- c(2e-1, 7e-1)

corr <- copula::p2P(rho[[1]], d = d)
walk(partition, ~{
  corr[.x, .x] <<- copula::p2P(rho[[2]], d = length(.x))
})

test_that("`simulate_dt` spearman correlation is in appr. CI for H2-Ext. Param", {
  expect_cor_not_rejected(
    simulate_dt(CuadrasAugeH2ExtMO3FParam(
      composition = composition, fraction = fraction, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
  expect_cor_not_rejected(
    simulate_dt(AlphaStableH2ExtMO3FParam(
      composition = composition, fraction = fraction, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
  expect_cor_not_rejected(
    simulate_dt(PoissonH2ExtMO3FParam(
      composition = composition, fraction = fraction, lambda = lambda, rho = rho), n_sim = n_sim),
    corr, level = level)
  expect_cor_not_rejected(
    simulate_dt(ExponentialH2ExtMO3FParam(
      composition = composition, fraction = fraction, lambda = lambda, rho = rho), n_sim = n_sim),
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
