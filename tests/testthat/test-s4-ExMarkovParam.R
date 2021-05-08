d <- 5L
lambda <- 8e-2
alpha <- 4e-1

bf <- ScaledBernsteinFunction(
  scale = lambda, original = AlphaStableBernsteinFunction(alpha = log2(2 - alpha)))
ex_qmatrix <- rmo::exQMatrix(bf, d = d)

test_that("`ExMarkovParam`-class is correctly initialized", {
  parm <- ExMarkovParam()
  expect_s4_class(parm, "ExMarkovParam")

  setDimension(parm) <- d
  setExQMatrix(parm) <- ex_qmatrix
  expect_true(isTRUE(validObject(parm)))
  expect_equal(getDimension(parm), d)
  expect_equal(getExQMatrix(parm), ex_qmatrix)

  expect_equal(parm, ExMarkovParam(ex_qmatrix))
})

test_that("`simulate_dt` works as expected for `ExMarkovParam`-class", {
  # HELPER START
  rfn <- function(n, ex_qmatrix) {
    qassert(n, "X1(0,)")
    assert_exqmatrix(ex_qmatrix)
    d <- nrow(ex_qmatrix) - 1L
    out <- matrix(nrow = n, ncol = d)
    for (k in 1:n) {
      state <- 0L
      time <- 0
      while (state != d) {
        wt <- rexp(1, rate = -ex_qmatrix[1L+state, 1L+state])
        time <- time + wt
        out[k, (1L+state):d] <- time
        state <- state +
          sample.int(n = d - state, size = 1L, replace = FALSE,
                      prob = ex_qmatrix[1L+state, (2L+state):(d+1L)])
      }
      perm <- sample.int(n = d, size = d, replace = FALSE)
      out[k, perm] <- out[k, perm]
    }
    if (isTRUE(nrow(out) <= 2L || ncol(out) <= 2L)) out <- as.vector(out)

    out
  }
  # HELPER END

  parm <- ExMarkovParam(ex_qmatrix)

  # n is 1, d is larger than 1
  set.seed(1623)
  x <- simulate_dt(parm, n_sim = 1L)
  expect_equal(length(x), d)

  set.seed(1623)
  y <- rfn(1L, ex_qmatrix)
  expect_equal(x, y)

  # n and d are larger than 1
  n <- 5e1L

  set.seed(1623)
  x <- simulate_dt(parm, n_sim = n)
  expect_equal(ncol(x), d)
  expect_equal(nrow(x), n)

  set.seed(1623)
  y <- rfn(n, ex_qmatrix)
  expect_equal(x, y)
})

test_that("`probability_distribution` works as expected for `ExMarkovParam`", {
  # HELPER START
  pfn <- function(t, ex_qmatrix) {
    qassert(t, "N+[0,)")
    assert_exqmatrix(ex_qmatrix)
    out <- purrr::map(t, ~{
        t(expm(.x * ex_qmatrix)[1L, , drop = FALSE])
      }) %>%
      purrr::reduce(cbind)
    if (isTRUE(nrow(out) <= 2L || ncol(out) <= 2L)) out <- as.vector(out)

    out
  }
  # HELPER END

  parm <- ExMarkovParam(ex_qmatrix)
  times <- seq(25e-2, 5L, by = 25e-2)

  # length of `times` is 1
  x <- probability_distribution(parm, times[[1]])
  expect_numeric(
    x, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, len = d+1L)
  expect_equal(sum(x), 1)
  expect_equal(x, pfn(times[[1]], ex_qmatrix))

  # length of `times` is larger than 1
  x <- probability_distribution(parm, times)
  expect_matrix(
    x, mode = "numeric", any.missing = FALSE, nrows = d+1L, ncols = length(times))
  expect_numeric(
    x, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE)
  expect_equal(apply(x, 2, sum), rep(1, ncol(x)))
  expect_equal(x, pfn(times, ex_qmatrix))

  # specify class explicitly
  expect_equal(x, probability_distribution(parm, times, method = "ExMarkovParam"))
})
