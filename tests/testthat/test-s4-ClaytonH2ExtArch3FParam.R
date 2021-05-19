partition <- list(1L:2L, 3L:6L, 7L:8L)
composition <- purrr::map_int(partition, length)
lambda <- 8e-2
rho <- c(0.25, 0.55)

d <- sum(composition)
nu <- purrr::map_dbl(rho, copula::iRho, copula = copula::claytonCopula())
tau <- purrr::map(nu, copula::claytonCopula) %>%
  purrr::map_dbl(copula::tau)

test_that("`ClaytonH2ExtArch3FParam`-class is correctly initialized", {
  parm <- ClaytonH2ExtArch3FParam()
  expect_s4_class(parm, "ClaytonH2ExtArch3FParam")

  setComposition(parm) <- composition
  setLambda(parm) <- lambda
  setNu(parm) <- nu
  expect_error(validObject(parm), NA)
  expect_equal(getDimension(parm), d)
  expect_equal(getComposition(parm), composition)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho, tolerance = testthat_tolerance()^0.5)
  expect_equal(getTau(parm), tau)

  expect_equal(parm, ClaytonH2ExtArch3FParam(composition, lambda, nu))
  expect_equal(parm, ClaytonH2ExtArch3FParam(composition, lambda, rho = rho))
  expect_equal(parm, ClaytonH2ExtArch3FParam(composition, lambda, tau = tau))
})

test_that("`ClaytonH2ExtArch3FParam`-class setters can be used in arbitrary order", { # nolint
  parm <- ClaytonH2ExtArch3FParam(composition, lambda, nu, fraction)

  parm2 <- ClaytonH2ExtArch3FParam()
  setComposition(parm2) <- composition
  setLambda(parm2) <- lambda
  setRho(parm2) <- rho
  expect_equal(parm, parm2)

  parm2 <- ClaytonH2ExtArch3FParam()
  setComposition(parm2) <- composition
  setLambda(parm2) <- lambda
  setTau(parm2) <- tau
  expect_equal(parm, parm2)

  parm2 <- ClaytonH2ExtArch3FParam()
  setComposition(parm2) <- composition
  setRho(parm2) <- rho
  setLambda(parm2) <- lambda
  expect_equal(parm, parm2)

  parm2 <- ClaytonH2ExtArch3FParam()
  setLambda(parm2) <- lambda
  setComposition(parm2) <- composition
  setRho(parm2) <- rho
  expect_equal(parm, parm2)

  parm2 <- ClaytonH2ExtArch3FParam()
  setRho(parm2) <- rho
  setComposition(parm2) <- composition
  setLambda(parm2) <- lambda
  expect_equal(parm, parm2)

  parm2 <- ClaytonH2ExtArch3FParam()
  setLambda(parm2) <- lambda
  setRho(parm2) <- rho
  setComposition(parm2) <- composition
  expect_equal(parm, parm2)

  parm2 <- ClaytonH2ExtArch3FParam()
  setRho(parm2) <- rho
  setLambda(parm2) <- lambda
  setComposition(parm2) <- composition
  expect_equal(parm, parm2)
})
