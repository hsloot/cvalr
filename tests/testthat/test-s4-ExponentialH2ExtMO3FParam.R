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
  expect_error(validObject(parm), NA)
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
  expect_equal(parm, ExponentialH2ExtMO3FParam(composition, lambda, rho = rho))
  expect_equal(parm, ExponentialH2ExtMO3FParam(composition, lambda, tau = tau))
  expect_equal(parm, ExponentialH2ExtMO3FParam(composition, lambda, alpha = alpha))
  expect_equal(as(parm, "H2ExMarkovParam"), H2ExMarkovParam(fraction, models))
  expect_equal(as(parm, "H2ExMOParam"), H2ExMOParam(fraction, models))
  expect_equal(as(parm, "H2ExtMOParam"), H2ExtMOParam(fraction, models))
})

test_that("`ExponentialH2ExtMO3FParam`-class setters can be used in arbitrary order", { # nolint
  parm <- ExponentialH2ExtMO3FParam(composition, lambda, nu, fraction)

  parm2 <- ExponentialH2ExtMO3FParam()
  setComposition(parm2) <- composition
  setLambda(parm2) <- lambda
  setRho(parm2) <- rho
  expect_equal(parm, parm2)

  parm2 <- ExponentialH2ExtMO3FParam()
  setComposition(parm2) <- composition
  setLambda(parm2) <- lambda
  setTau(parm2) <- tau
  expect_equal(parm, parm2)

  parm2 <- ExponentialH2ExtMO3FParam()
  setComposition(parm2) <- composition
  setRho(parm2) <- rho
  setLambda(parm2) <- lambda
  expect_equal(parm, parm2)

  parm2 <- ExponentialH2ExtMO3FParam()
  setLambda(parm2) <- lambda
  setComposition(parm2) <- composition
  setRho(parm2) <- rho
  expect_equal(parm, parm2)

  parm2 <- ExponentialH2ExtMO3FParam()
  setRho(parm2) <- rho
  setComposition(parm2) <- composition
  setLambda(parm2) <- lambda
  expect_equal(parm, parm2)

  parm2 <- ExponentialH2ExtMO3FParam()
  setLambda(parm2) <- lambda
  setRho(parm2) <- rho
  setComposition(parm2) <- composition
  expect_equal(parm, parm2)

  parm2 <- ExponentialH2ExtMO3FParam()
  setRho(parm2) <- rho
  setLambda(parm2) <- lambda
  setComposition(parm2) <- composition
  expect_equal(parm, parm2)
})
