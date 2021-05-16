partition <- list(1L:2L, 3L:6L, 7L:8L)
composition <- purrr::map_int(partition, length)
lambda <- 8e-2
rho <- c(0.25, 0.55)

d <- sum(composition)
nu <- purrr::map_dbl(rho, copula::iRho, copula = copula::gumbelCopula())
tau <- purrr::map(nu, copula::gumbelCopula) %>%
  purrr::map_dbl(copula::tau)

test_that("`GumbelH2ExtArch3FParam`-class is correctly initialized", {
  parm <- GumbelH2ExtArch3FParam()
  expect_s4_class(parm, "GumbelH2ExtArch3FParam")

  setComposition(parm) <- composition
  setLambda(parm) <- lambda
  setNu(parm) <- nu
  expect_true(validObject(parm))
  expect_equal(getDimension(parm), d)
  expect_equal(getComposition(parm), composition)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho, tolerance = testthat_tolerance()^0.5)
  expect_equal(getTau(parm), tau)

  expect_equal(parm, GumbelH2ExtArch3FParam(composition, lambda, nu))
  expect_equal(parm, GumbelH2ExtArch3FParam(composition, lambda, rho = rho))
  expect_equal(parm, GumbelH2ExtArch3FParam(composition, lambda, tau = tau))
})
