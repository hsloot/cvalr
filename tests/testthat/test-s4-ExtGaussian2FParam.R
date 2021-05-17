d <- 5L
lambda <- 8e-2
nu <- 4e-1

rho <- (6 / pi) * asin(nu / 2)
tau <- (2 / pi) * asin(nu)


test_that("`ExtGaussian2FParam`-class is correctly initialized", {
  parm <- ExtGaussian2FParam()
  expect_s4_class(parm, "ExtGaussian2FParam")

  setDimension(parm) <- d
  setLambda(parm) <- lambda
  setNu(parm) <- nu
  expect_error(validObject(parm), NA)
  expect_equal(getDimension(parm), d)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)

  expect_equal(parm, ExtGaussian2FParam(d, lambda, nu))
  expect_equal(parm, ExtGaussian2FParam(d, lambda, rho = rho))
  expect_equal(parm, ExtGaussian2FParam(d, lambda, tau = tau))
})

test_that("`ExtGaussian2FParam`-class setters can be used in arbitrary order", { # nolint
  parm <- ExtGaussian2FParam(d, lambda, nu)

  parm2 <- ExtGaussian2FParam()
  setDimension(parm2) <- d
  setLambda(parm2) <- lambda
  setRho(parm2) <- rho
  expect_equal(parm, parm2)

  parm2 <- ExtGaussian2FParam()
  setDimension(parm2) <- d
  setLambda(parm2) <- lambda
  setTau(parm2) <- tau
  expect_equal(parm, parm2)

  parm2 <- ExtGaussian2FParam()
  setDimension(parm2) <- d
  setNu(parm2) <- nu
  setLambda(parm2) <- lambda
  expect_equal(parm, parm2)

  parm2 <- ExtGaussian2FParam()
  setLambda(parm2) <- lambda
  setDimension(parm2) <- d
  setNu(parm2) <- nu
  expect_equal(parm, parm2)

  parm2 <- ExtGaussian2FParam()
  setNu(parm2) <- nu
  setDimension(parm2) <- d
  setLambda(parm2) <- lambda
  expect_equal(parm, parm2)

  parm2 <- ExtGaussian2FParam()
  setLambda(parm2) <- lambda
  setNu(parm2) <- nu
  setDimension(parm2) <- d
  expect_equal(parm, parm2)

  parm2 <- ExtGaussian2FParam()
  setNu(parm2) <- nu
  setLambda(parm2) <- lambda
  setDimension(parm2) <- d
  expect_equal(parm, parm2)
})

test_that("`simulate_dt` works as expected for `ExtGaussian2FParam`", {
  # HELPER START
  rfn <- function(parm, n) {
    qassert(n, "X1(0,)")
    copula::normalCopula(param = getNu(parm), dim = getDimension(parm), dispstr = "ex") %>%
      copula::rCopula(n, .) %>%
      qexp(rate = getLambda(parm), lower.tail = FALSE)
  }
  # HELPER END

  parm <- ExtGaussian2FParam(d, lambda, nu)

  # n is 1, d is larger than 1
  set.seed(1623)
  x <- simulate_dt(parm, n_sim = 1L)
  expect_numeric(
    x, lower = 0, finite = TRUE, any.missing = FALSE, len = d)

  set.seed(1623)
  y <- rfn(parm, 1L)
  expect_equal(x, y)

  # n and d are larger than 1
  n <- 5e1L

  set.seed(1623)
  x <- simulate_dt(parm, n_sim = n)
  expect_matrix(
    x, mode = "numeric", any.missing = FALSE, nrows = n, ncols = d)
  expect_numeric(
    x, lower = 0, finite = TRUE)

  set.seed(1623)
  y <- rfn(parm, n)
  expect_equal(x, y)
})


test_that("`probability_distribution` works as expected for `ExtGaussian2FParam`", {
  # HELPER START
  pfn <- function(parm, times) {
    qassert(times, "N+[0,)")

    d <- getDimension(parm)
    lambda <- getLambda(parm)
    nu <- getNu(parm)
    out <- matrix(nrow = d+1L, ncol = length(times))
    for (j in seq_along(times)) {
      if (0 == times[j]) {
        out[, j] <- c(1, rep(0, d))
      } else if (Inf == times[j]) {
        out[, j] <- c(rep(0, d), 1)
      } else {
        for (k in 0L:d) {
          int_fn <- function(x) {
            p <- pnorm((qnorm(pexp(times[j], rate = lambda)) - sqrt(nu) * x) / (sqrt(1 - nu)))
            v_multiply_binomial_coefficient(p^k * (1 - p) ^ (d-k) * dnorm(x), d, k)
          }
          int_res <- integrate(
            int_fn, lower = -Inf, upper = Inf, rel.tol = .Machine$double.eps^0.5)
          out[k+1L, j] <- int_res$value
        }
      }
    }

    out
  }
  # HELPER END

  parm <- ExtGaussian2FParam(d, lambda, nu)
  times <- seq(25e-2, 5L, by = 25e-2)

  # length of `times` is 1
  x <- probability_distribution(parm, times[[1]])
  expect_numeric(
    x, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, len = d+1L)
  expect_equal(sum(x), 1)
  expect_equal(x, pfn(parm, times[[1]]))

  # length of `times` is larger than 1
  x <- probability_distribution(parm, times)
  expect_matrix(
    x, mode = "numeric", any.missing = FALSE, nrows = d+1L, ncols = length(times))
  expect_numeric(
    x, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE)
  expect_equal(apply(x, 2, sum), rep(1, ncol(x)))
  expect_equal(x, pfn(parm, times))

  # specify class explicitly
  expect_equal(x, probability_distribution(parm, times, method = "ExtGaussian2FParam"))
})


test_that("`expected_pcds_loss` works as expected for `ExtGaussian2FParam", {
  # HELPER START
  epcdslfn <- function(parm, times, recovery_rate) {
    qassert(times, "N+[0,)")
    qassert(recovery_rate, "N1[0,1]")

    (1 - recovery_rate) * pexp(times, rate = getLambda(parm))
  }
  # HELPER END

  parm <- ExtGaussian2FParam(dim = d, lambda = lambda, nu = nu)
  times <- seq(25e-2, 5L, by = 25e-2)
  recovery_rate <- 0.4

  x <- expected_pcds_loss(parm, times, recovery_rate = recovery_rate)
  expect_numeric(x, any.missing = FALSE, lower = 0, upper = 1,
    len = length(times), sorted = TRUE)
  expect_equal(x, epcdslfn(parm, times, recovery_rate))
})


test_that("`expected_cdo_loss` works as expected for `ExtGaussian2FParam`", {
  # HELPER START
  ecdolfn <- function(parm, times, recovery_rate, lower, upper) {
    qassert(times, "N+[0,)")
    qassert(recovery_rate, "N1[0,1]")
    qassert(lower, "N1[0,1]")
    qassert(upper, "N1[0,1]")
    lambda <- getLambda(parm)
    nu <- getNu(parm)

    cop <- copula::normalCopula(param = -sqrt(1 - nu), dim = 2L, dispstr = "ex")
    u_left <- cbind(1 - pmin(lower / (1 - recovery_rate), 1), pexp(times, rate = lambda))
    u_right <- cbind(1 - pmin(upper / (1 - recovery_rate), 1), pexp(times, rate = lambda))
    (1 - recovery_rate) * (copula::pCopula(u_left, cop) - copula::pCopula(u_right, cop))
  }
  # HELPER END

  parm <- ExtGaussian2FParam(dim = d, lambda = lambda, nu = nu)
  times <- seq(25e-2, 5L, by = 25e-2)
  recovery_rate <- 0.4
  lower <- 0.1
  upper <- 0.2

  x <- expected_cdo_loss(parm, times, recovery_rate, lower, upper)
  expect_numeric(x, any.missing = FALSE, lower = 0, upper = 1,
    len = length(times), sorted = TRUE)
  expect_equal(x, ecdolfn(parm, times, recovery_rate, lower, upper))
})
