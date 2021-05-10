partition <- list(1L:2L, 3L:6L, 7L:8L)
lambda <- 8e-2
nu <- c(0.3, 0.8)

d <- length(unlist(partition))
rho <- (6 / pi) * asin(nu / 2)
tau <- (2 / pi) * asin(nu)

test_that("`H2ExtGaussian3FParam`-class is correctly initialized", {
  parm <- H2ExtGaussian3FParam()
  expect_s4_class(parm, "H2ExtGaussian3FParam")

  setPartition(parm) <- partition
  setLambda(parm) <- lambda
  setNu(parm) <- nu
  expect_true(validObject(parm))
  expect_equal(getDimension(parm), d)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)

  expect_equal(parm, H2ExtGaussian3FParam(partition, lambda, nu))
  expect_equal(parm, H2ExtGaussian3FParam(partition, lambda, rho = rho))
  expect_equal(parm, H2ExtGaussian3FParam(partition, lambda, tau = tau))
})
