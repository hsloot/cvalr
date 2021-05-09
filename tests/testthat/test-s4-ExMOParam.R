d <- 5L
lambda <- 8e-2
alpha <- 4e-1

bf <- ScaledBernsteinFunction(
  scale = lambda, original = AlphaStableBernsteinFunction(alpha = log2(2 - alpha)))
ex_qmatrix <- rmo::exQMatrix(bf, d)
ex_intensities <- rmo::exIntensities(bf, d = d)

test_that("`ExMOParam`-class is correctly initialized", {
  parm <- ExMOParam()
  expect_s4_class(parm, "ExMOParam")

  setExIntensities(parm) <- ex_intensities
  expect_true(isTRUE(validObject(parm)))
  expect_equal(getDimension(parm), d)
  expect_equal(getExQMatrix(parm), ex_qmatrix)
  expect_equal(getExIntensities(parm), ex_intensities)

  expect_equal(parm, ExMOParam(ex_intensities))
  expect_equal(as(parm, "ExMarkovParam"), ExMarkovParam(ex_qmatrix))
})

test_that("`simulate_dt` works as expected for `ExMOParam`", {
  # HELPER START
  rfn <- function(n, ex_intensities) {
    qassert(n, "X1(0,)")
    assert_numeric(
      ex_intensities, lower = 0, upper = Inf, finite = TRUE, any.missing = FALSE)
    d <- length(ex_intensities)
    out <- rmo::rexmo_markovian(n, d, ex_intensities)
    if (isTRUE(nrow(out) <= 2L || ncol(out) <= 2L)) out <- as.vector(out)

    out
  }
  # HELPER END

  parm <- ExMOParam(ex_intensities)

  # n is 1, d is larger than 1
  set.seed(1623)
  x <- simulate_dt(parm, n_sim = 1L)
  expect_numeric(
    x, lower = 0, finite = TRUE, any.missing = FALSE, len = d)

  set.seed(1623)
  y <- rfn(1L, ex_intensities)
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
  y <- rfn(n, ex_intensities)
  expect_equal(x, y)
})
