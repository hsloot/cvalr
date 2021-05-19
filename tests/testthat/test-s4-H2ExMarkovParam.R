partition <- list(1L:2L, 3L:6L, 7L:8L)
composition <- purrr::map_int(partition, length)
lambda <- 8e-2
alpha <- c(0.3, 0.8)

d <- sum(composition)
rho <- 3 * alpha / (4 - alpha)
tau <- alpha / (2 - alpha)

fraction <- (2 * alpha[[1]] + 1 - alpha[[2]]) / 2
models <- purrr::map(composition, ~{
  as(
    AlphaStableExtMO2FParam(.x, lambda, alpha = (alpha[[2]] - alpha[[1]]) / (1 - fraction)),
    "ExMarkovParam")
})
models <- c(
  list(as(AlphaStableExtMO2FParam(d, lambda, alpha = alpha[[1]] / fraction), "ExMarkovParam")),
  models)

test_that("`H2ExMarkovParam`-class is correctly initialized", {
  parm <- H2ExMarkovParam()
  expect_s4_class(parm, "H2ExMarkovParam")

  setFraction(parm) <- fraction
  setModels(parm) <- models
  expect_error(validObject(parm), NA)
  expect_equal(getDimension(parm), d)
  expect_equal(getComposition(parm), composition)
  expect_equal(getPartition(parm), partition)
  expect_equal(getFraction(parm), fraction)
  expect_equal(getModels(parm), models)
  expect_equal(getGlobalModel(parm), models[[1L]])
  expect_equal(getPartitionModels(parm), models[-1L])

  expect_equal(parm, H2ExMarkovParam(fraction, models))
})

test_that("`simulate_dt` works as expected for `H2ExMarkovParam`-class", {
  # HELPER START
  rfn <- function(parm, n) {
    qassert(n, "X1[1,)")
    d <- getDimension(parm)
    fraction <- getFraction(parm)
    parm_global <- getGlobalModel(parm)
    parm_partition <- getPartitionModels(parm)
    tmp_global <- simulate_dt(parm_global, n_sim = n)
    tmp_composition <- purrr::map(parm_partition, simulate_dt, n_sim = n) %>%
        purrr::reduce(cbind) %>%
        `dimnames<-`(NULL)

    pmin(1 / fraction * tmp_global, 1 / (1 - fraction) * tmp_composition)
  }
  # HELPER END

  parm <- H2ExMarkovParam(fraction, models)

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
