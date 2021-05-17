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
    "ExtMOParam")
})
models <- c(
  list(as(AlphaStableExtMO2FParam(d, lambda, alpha = alpha[[1]] / fraction), "ExtMOParam")),
  models)

test_that("`H2ExtMOParam`-class is correctly initialized", {
  parm <- H2ExtMOParam()
  expect_s4_class(parm, "H2ExtMOParam")

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

  expect_equal(parm, H2ExtMOParam(fraction, models))
  expect_equal(as(parm, "H2ExMarkovParam"), H2ExMarkovParam(fraction, models))
  expect_equal(as(parm, "H2ExMOParam"), H2ExMOParam(fraction, models))
})
