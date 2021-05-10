d <- 5L
lambda <- 8e-2
rho <- 4e-1


n <- 5e1L

test_that("`simulate_dt` is working as expected for `ExtArch2FParam`", {
  # HELPER START
  rfn <- function(n, parm) {
    qassert(n, "X1(0,)")
    assert_class(parm, "ExtArch2FParam")
    assert_true(validObject(parm))

    out <- qexp(
      copula::rCopula(n, parm@copula), rate = lambda, lower.tail = !parm@survival)
      if (isTRUE(nrow(out) == 1L || ncol(out) == 1L)) out <- as.vector(out)

      out
  }
  # HELPER END
  parm <- FrankExtArch2FParam(dim = d, lambda = lambda, rho = rho)

  set.seed(1623)
  x <- simulate_dt(parm, n_sim = n)
  expect_matrix(
    x, mode = "numeric", any.missing = FALSE, nrows = n, ncols = d)
  expect_numeric(x, lower = 0, finite = TRUE)

  set.seed(1623)
  expect_equal(x, rfn(n, parm))
})


times <- seq(0, 5, by = 0.25)
recovery_rate <- 4e-1

test_that("`expected_pcds_loss` works as expected for `ExtArch2FParam", {
  parm <- FrankExtArch2FParam(dim = d, lambda = lambda, rho = rho)

  x <- expected_pcds_loss(parm, times, recovery_rate = recovery_rate)
  expect_numeric(
    x, any.missing = FALSE, lower = 0, upper = 1, len = length(times), sorted = TRUE)
  expect_equal(x, (1 - recovery_rate) * pexp(times, rate = lambda))
})
