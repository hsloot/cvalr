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

test_that("`expected_pcds_equation` works as expected for `ExtArch2FParam", {
  times <- seq(25e-2, 5, by = 25e-2)
  discount_factors <- rep(1, length(times))
  recovery_rate <- 4e-1
  coupon <- 1e-1
  upfront <- -1e-2
  parm <- FrankExtArch2FParam(dim = d, lambda = lambda, rho = rho)

  x <- expected_pcds_equation(parm, times, discount_factors, recovery_rate, coupon, upfront)
  expect_numeric(x, finite = TRUE, any.missing = FALSE, len = 1L)
  y <- test__expected_pcds_equation__default(
    parm, times, discount_factors, recovery_rate, coupon, upfront)
  expect_equal(x, y)
})
