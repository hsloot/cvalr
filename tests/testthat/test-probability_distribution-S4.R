times <- seq(0, 5, by = 0.25)

dim <- 50L
lambda <- 8e-2
rho <- 4e-1

seed <- 1623
n_sim <- 1e2


test_that("`probability_distribution` works as expected for `CalibrationParam`", {
  parm <- AlphaStableExtMO2FParam(dim = dim, lambda = lambda, rho = rho)

  x <- probability_distribution(parm, times,
        method = "CalibrationParam", seed = 1623, sim_args = list(n_sim = n_sim))
  expect_matrix(x,
    mode = "numeric", any.missing = FALSE, nrows = dim+1, ncols = length(times))
  expect_numeric(x,
    any.missing = FALSE, lower = 0, upper = 1)
  expect_equal(as.vector(apply(x, 2, sum)), rep(1, ncol(x)))

  set.seed(seed)
  expect_equal(x,
    probability_distribution(parm, times,
      method = "CalibrationParam", sim_args = list(n_sim = n_sim)))

  pd_naive <- function(object, times, n_sim) {
    x <- simulate_adcp(object, times, n_sim = n_sim)
    if (!is.matrix(x)) {
      x <- as.matrix(x, ncol = 1L)
    }
    out <- t(sapply((0:object@dim) / object@dim, function(k) {
      apply(x, 2, function(.x) mean(.x == k))
    }))

    if (1L == nrow(out) || 1L == ncol(out)) {
      out <- as.vector(out)
    }

    out
  }

  set.seed(seed)
  expect_equal(x, pd_naive(parm, times, n_sim))
})


test_that("`probability_distribution` works as expected for `ExtGaussian2FParam`", { # nolint
  parm <- ExtGaussian2FParam(dim = dim, lambda = lambda, rho = rho)

  x <- probability_distribution(parm, times)
  expect_matrix(x,
    mode = "numeric", any.missing = FALSE, nrows = dim+1, ncols = length(times))
  expect_numeric(x,
    any.missing = FALSE, lower = 0, upper = 1)
  expect_equal(as.vector(apply(x, 2, sum)), rep(1, ncol(x)))

  pd_naive <- function(dim, lambda, nu, times) {
    out <- matrix(nrow = dim+1, ncol = length(times))
    times <- qnorm(pexp(times, rate = lambda))
    for (k in 0:dim) {
      for (j in seq_along(times)) {
        if (-Inf == times[j] && 0 == k) {
          out[1+k, j] <- 1
        } else if (-Inf == times[j] && 0 < k) {
          out[1+k, j] <- 0
        } else {
          out[1+k, j] <- integrate(
            function(x) {
              ldp <- pnorm(
                (times[j] - sqrt(nu) * x) / sqrt(1 - nu),
                log.p = TRUE, lower.tail = TRUE)
              lsp <- pnorm(
                (times[j] - sqrt(nu) * x) / sqrt(1 - nu),
                log.p = TRUE, lower.tail = FALSE)
              sapply(
                exp(k * ldp + (dim - k) * lsp) * dnorm(x),
                function(v) {
                  multiply_binomial_coefficient(v, dim, k)
                })
            }, lower = -Inf, upper = Inf, rel.tol = .Machine$double.eps^0.5)$value
        }
      }
    }

    out
  }

  expect_equal(x, pd_naive(parm@dim, parm@lambda, parm@nu, times))
})
