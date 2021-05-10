partition <- list(1L:2L, 3L:6L, 7L:8L)
lambda <- 8e-2
tau <- c(0.25, 0.55)

d <- length(unlist(partition))
nu <- purrr::map_dbl(tau, copula::iTau, copula = copula::joeCopula())

test_that("`JoeH2ExtArch3FParam`-class is correctly initialized", {
  parm <- JoeH2ExtArch3FParam()
  expect_s4_class(parm, "JoeH2ExtArch3FParam")

  setPartition(parm) <- partition
  setLambda(parm) <- lambda
  setNu(parm) <- nu
  expect_true(validObject(parm))
  expect_equal(getDimension(parm), d)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getTau(parm), tau)

  expect_equal(parm, JoeH2ExtArch3FParam(partition, lambda, nu))
  expect_equal(parm, JoeH2ExtArch3FParam(partition, lambda, tau = tau))
})
