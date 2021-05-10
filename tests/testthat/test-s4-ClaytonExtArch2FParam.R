d <- 5L
lambda <- 0.08
rho <- 0.45

nu <- copula::iRho(copula::claytonCopula(), rho)
cop <- copula::claytonCopula(param = nu)
tau <- copula::tau(cop)

test_that("`ClaytonExtArch2FParam`-class is correctly initialized", {
  parm <- ClaytonExtArch2FParam()
  expect_s4_class(parm, "ClaytonExtArch2FParam")

  parm@survival <- FALSE
  parm@copula <- copula::archmCopula(family = "Clayton")
  setDimension(parm) <- d
  setLambda(parm) <- lambda
  setNu(parm) <- nu
  expect_true(validObject(parm))
  expect_equal(getDimension(parm), d)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho, tolerance = .Machine$double.eps^0.25)
  expect_equal(getTau(parm), tau)

  expect_equal(parm, ClaytonExtArch2FParam(d, lambda, nu))
  expect_equal(parm, ClaytonExtArch2FParam(d, lambda, rho = rho))
  expect_equal(parm, ClaytonExtArch2FParam(d, lambda, tau = tau))
})
