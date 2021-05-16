partition <- list(1L:2L, 3L:6L, 7L:8L)
composition <- purrr::map_int(partition, length)
lambda <- 8e-2
alpha <- c(0.3, 0.8)

d <- sum(composition)
rho <- 3 * alpha / (4 - alpha)
tau <- alpha / (2 - alpha)

fraction <- (2 * alpha[[1]] + 1 - alpha[[2]]) / 2
models <- purrr::map(composition, ~{
  ExponentialExtMO2FParam(.x, lambda, alpha = (alpha[[2]] - alpha[[1]]) / (1 - fraction))
})
models <- c(
  list(ExponentialExtMO2FParam(d, lambda, alpha = alpha[[1]] / fraction)),
  models)
nu <- purrr::map_dbl(models[1:2], getNu)

test_that("`ExponentialH2ExtMO3FParam`-class is correctly initialized", {
  parm <- ExponentialH2ExtMO3FParam()
  expect_s4_class(parm, "ExponentialH2ExtMO3FParam")

  setComposition(parm) <- composition
  setFraction(parm) <- fraction
  setLambda(parm) <- lambda
  setNu(parm) <- nu
  expect_true(validObject(parm))
  expect_equal(getDimension(parm), d)
  expect_equal(getComposition(parm), composition)
  expect_equal(getPartition(parm), partition)
  expect_equal(getModels(parm), models)
  expect_equal(getFraction(parm), fraction)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)

  expect_equal(parm, ExponentialH2ExtMO3FParam(composition, lambda, nu, fraction))
  expect_equal(parm, ExponentialH2ExtMO3FParam(composition, lambda, fraction = fraction, rho = rho))
  expect_equal(parm, ExponentialH2ExtMO3FParam(composition, lambda, fraction = fraction, tau = tau))
  # TODO: Replace after implementation of `coerce`
  expect_equal(as(parm, "H2ExMarkovParam"), H2ExMarkovParam(models, fraction))
  expect_equal(as(parm, "H2ExMOParam"), H2ExMOParam(models, fraction))
  expect_equal(as(parm, "H2ExtMOParam"), H2ExtMOParam(models, fraction))
  # nolint start
  # expect_equal(as(parm, "H2ExMarkovParam"), H2ExMarkovParam(purrr::map(models, as, Class = "ExMarkovParam"), fraction))
  # expect_equal(as(parm, "H2ExMOParam"), H2ExMOParam(purrr::map(models, as, Class = "ExMOParam"), fraction))
  # expect_equal(as(parm, "H2ExtMOParam"), H2ExtMOParam(purrr::map(models, as, Class = "ExtMOParam"), fraction))
  # nolint end
})
