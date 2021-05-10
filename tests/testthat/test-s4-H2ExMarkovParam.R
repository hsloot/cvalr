partition <- list(1L:2L, 3L:6L, 7L:8L)
lambda <- 8e-2
alpha <- c(0.3, 0.8)

d <- length(unlist(partition))
rho <- 3 * alpha / (4 - alpha)
tau <- alpha / (2 - alpha)

fraction <- (2 * alpha[[1]] + 1 - alpha[[2]]) / 2
models <- purrr::map(partition, ~{
  .d <- length(.x)
  as(
    AlphaStableExtMO2FParam(.d, lambda, alpha = (alpha[[2]] - alpha[[1]]) / (1 - fraction)),
    "ExMarkovParam")
})
models <- c(
  list(as(AlphaStableExtMO2FParam(d, lambda, alpha = alpha[[1]] / fraction), "ExMarkovParam")),
  models)

test_that("`H2ExMarkovParam`-class is correctly initialized", {
  parm <- H2ExMarkovParam()
  expect_s4_class(parm, "H2ExMarkovParam")

  setModels(parm) <- models
  setFraction(parm) <- fraction
  expect_true(validObject(parm))
  expect_equal(getDimension(parm), d)
  expect_equal(getPartition(parm), partition)
  expect_equal(getModels(parm), models)
  expect_equal(getFraction(parm), fraction)

  expect_equal(parm, H2ExMarkovParam(models, fraction))
})
