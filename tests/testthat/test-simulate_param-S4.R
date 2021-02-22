times <- seq(0, 5, by = 0.25)

dim = 50L
lambda <- 8e-2
rho <- 4e-1

seed <- 1623
n_sim <- 1e2

dt2dct <- function(x, times) {
  out <- matrix(nrow = nrow(x), ncol = length(times))
  for (k in 1:nrow(out)) {
    for (j in 1:ncol(out)) {
      out[k, j] <- mean(x[k, ] <= times[j])
    }
  }

  out
}

test_that("`simulate_param` is working correctly for `ExMarkovParam`", {
  parm <- AlphaStableExtMO2FParam(dim = dim, lambda = lambda, rho = rho)

  set.seed(seed)
  x <- simulate_param(as(parm, "ExMarkovParam"), times, n_sim = n_sim)
  expect_matrix(x,
    mode = "numeric", any.missing = FALSE, nrows = n_sim, ncols = length(times))
  expect_integerish(x * parm@dim,
                any.missing = FALSE,
                lower = 0, upper = parm@dim)
  expect_matrix(x,
                mode = "numeric", any.missing = FALSE,
                nrows = n_sim, ncols = length(times))
  expect_integerish(t(apply(x * parm@dim, 1, diff)),
    any.missing = FALSE, lower = 0, upper = parm@dim)

  set.seed(seed)
  expect_equal(x,
    simulate_param(parm, times, method = "ExMarkovParam", n_sim = n_sim))

  sample_naive <- function(n, qmatrix, times) {
    dim <- nrow(qmatrix) - 1L
    tmp <- matrix(nrow = n_sim, ncol = dim)
    for (k in 1:n_sim) {
      state <- 0
      time <- 0
      while (state != dim) {
        wt <- rexp(1, rate = -qmatrix[1+state, 1+state])
        time <- time + wt
        tmp[k, (1+state):dim] <- time
        state <- state +
          sample.int(n = dim-state, size = 1, replace = FALSE,
            prob = qmatrix[1+state, (2+state):(dim+1)])
      }
    }
    out <- dt2dct(tmp, times)

    if (1L == nrow(out) || 1L == ncol(out))
      out <- as.vector(out)

    out
  }

  set.seed(seed)
  expect_equal(x, sample_naive(n_sim, parm@qmatrix, times))
})


test_that("`simulate_param` is working correctly for `ExMOParam`", {
  parm <- AlphaStableExtMO2FParam(dim = dim, lambda = lambda, rho = rho)

  set.seed(seed)
  x <- simulate_param(parm, times, n_sim = n_sim)
  expect_matrix(x,
                mode = "numeric", any.missing = FALSE, nrows = n_sim, ncols = length(times))
  expect_integerish(x * parm@dim,
                    any.missing = FALSE,
                    lower = 0, upper = parm@dim)
  expect_matrix(x,
                mode = "numeric", any.missing = FALSE,
                nrows = n_sim, ncols = length(times))
  expect_integerish(t(apply(x * parm@dim, 1, diff)),
                    any.missing = FALSE, lower = 0, upper = parm@dim)

  sample_naive <- function(n, ex_intensities, times) {
    dim <- length(ex_intensities)
    tmp <- rmo::rexmo_markovian(n_sim, dim, ex_intensities)
    out <- dt2dct(tmp, times)

    if (1L == nrow(out) || 1L == ncol(out))
      out <- as.vector(out)

    out
  }

  set.seed(seed)
  expect_equal(x, sample_naive(n_sim, parm@ex_intensities, times))
})


test_that("`simulate_param` is working correctly for `ExtGaussian2FParam`", {
  parm <- ExtGaussian2FParam(dim = dim, lambda = lambda, rho = rho)

  set.seed(seed)
  x <- simulate_param(parm, times, n_sim = n_sim)
  expect_matrix(x,
                mode = "numeric", any.missing = FALSE, nrows = n_sim, ncols = length(times))
  expect_integerish(x * parm@dim,
                    any.missing = FALSE,
                    lower = 0, upper = parm@dim)
  expect_matrix(x,
                mode = "numeric", any.missing = FALSE,
                nrows = n_sim, ncols = length(times))
  expect_integerish(t(apply(x * parm@dim, 1, diff)),
                    any.missing = FALSE, lower = 0, upper = parm@dim)

  sample_naive <- function(n, dim, lambda, nu, times) {
    tmp <- stats::qexp(
      copula::rCopula(n_sim, normalCopula(nu, dim = dim, dispstr = "ex")),
      rate = lambda, lower.tail = FALSE)
    out <- dt2dct(tmp, times)

    if (1L == nrow(out) || 1L == ncol(out))
      out <- as.vector(out)

    out
  }

  set.seed(seed)
  expect_equal(x, sample_naive(n_sim, parm@dim, parm@lambda, parm@nu, times))
})


test_that("`simulate_param` is working correctly for `FrankExtArch2FParam`", {
  parm <- FrankExtArch2FParam(dim = dim, lambda = lambda, rho = rho)

  set.seed(seed)
  x <- simulate_param(parm, times, n_sim = n_sim)
  expect_matrix(x,
                mode = "numeric", any.missing = FALSE, nrows = n_sim, ncols = length(times))
  expect_integerish(x * parm@dim,
                    any.missing = FALSE,
                    lower = 0, upper = parm@dim)
  expect_matrix(x,
                mode = "numeric", any.missing = FALSE,
                nrows = n_sim, ncols = length(times))
  expect_integerish(t(apply(x * parm@dim, 1, diff)),
                    any.missing = FALSE, lower = 0, upper = parm@dim)

  sample_naive <- function(n, dim, lambda, nu, times) {
    tmp <- stats::qexp(
      copula::rCopula(n_sim, frankCopula(nu, dim = dim)),
      rate = lambda, lower.tail = FALSE)
    out <- dt2dct(tmp, times)

    if (1L == nrow(out) || 1L == ncol(out))
      out <- as.vector(out)

    out
  }

  set.seed(seed)
  expect_equal(x, sample_naive(n_sim, parm@dim, parm@lambda, parm@nu, times))
})
