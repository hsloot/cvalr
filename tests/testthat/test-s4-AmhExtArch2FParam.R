d <- 5L
lambda <- 0.08
rho <- 0.45

nu <- copula::iRho(copula::amhCopula(), rho)
cop <- copula::amhCopula(param = nu)
tau <- copula::tau(cop)

test_that("`AmhExtArch2FParam`-class is correctly initialized", {
  parm <- AmhExtArch2FParam()
  expect_s4_class(parm, "AmhExtArch2FParam")

  parm@survival <- TRUE
  parm@copula <- copula::archmCopula(family = "Amh")
  setDimension(parm) <- d
  setLambda(parm) <- lambda
  setNu(parm) <- nu
  expect_true(validObject(parm))
  expect_equal(getDimension(parm), d)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho, tolerance = .Machine$double.eps^0.25)
  expect_equal(getTau(parm), tau)

  expect_equal(parm, AmhExtArch2FParam(d, lambda, nu))
  expect_equal(parm, AmhExtArch2FParam(d, lambda, rho = rho))
  expect_equal(parm, AmhExtArch2FParam(d, lambda, tau = tau))
})
