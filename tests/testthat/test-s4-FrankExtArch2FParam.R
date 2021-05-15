d <- 4L
lambda <- 0.08
rho <- 0.45

nu <- copula::iRho(copula::frankCopula(), rho)
cop <- copula::frankCopula(param = nu)
tau <- copula::tau(cop)

test_that("`FrankExtArch2FParam`-class is correctly initialized", {
  parm <- FrankExtArch2FParam()
  expect_s4_class(parm, "FrankExtArch2FParam")

  parm@survival <- TRUE
  parm@copula <- copula::archmCopula(family = "Frank")
  setDimension(parm) <- d
  setLambda(parm) <- lambda
  setNu(parm) <- nu
  expect_true(validObject(parm))
  expect_equal(getDimension(parm), d)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho, tolerance = .Machine$double.eps^0.25)
  expect_equal(getTau(parm), tau)

  expect_equal(parm, FrankExtArch2FParam(d, lambda, nu))
  expect_equal(parm, FrankExtArch2FParam(d, lambda, rho = rho))
  expect_equal(parm, FrankExtArch2FParam(d, lambda, tau = tau))
})

test_that("`FrankExtArch2FParam`-class setters can be used in arbitrary order", { # nolint
  parm <- FrankExtArch2FParam(d, lambda, nu)

  parm2 <- FrankExtArch2FParam()
  setDimension(parm2) <- d
  setLambda(parm2) <- lambda
  setRho(parm2) <- rho
  expect_equal(parm, parm2)

  parm2 <- FrankExtArch2FParam()
  setDimension(parm2) <- d
  setLambda(parm2) <- lambda
  setTau(parm2) <- tau
  expect_equal(parm, parm2)

  parm2 <- FrankExtArch2FParam()
  setDimension(parm2) <- d
  setNu(parm2) <- nu
  setLambda(parm2) <- lambda
  expect_equal(parm, parm2)

  parm2 <- FrankExtArch2FParam()
  setLambda(parm2) <- lambda
  setDimension(parm2) <- d
  setNu(parm2) <- nu
  expect_equal(parm, parm2)

  parm2 <- FrankExtArch2FParam()
  setNu(parm2) <- nu
  setDimension(parm2) <- d
  setLambda(parm2) <- lambda
  expect_equal(parm, parm2)

  parm2 <- FrankExtArch2FParam()
  setLambda(parm2) <- lambda
  setNu(parm2) <- nu
  setDimension(parm2) <- d
  expect_equal(parm, parm2)

  parm2 <- FrankExtArch2FParam()
  setNu(parm2) <- nu
  setLambda(parm2) <- lambda
  setDimension(parm2) <- d
  expect_equal(parm, parm2)
})
