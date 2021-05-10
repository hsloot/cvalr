d <- 5L
lambda <- 0.08
tau <- 0.25

nu <- copula::iTau(copula::joeCopula(), tau)
cop <- copula::joeCopula(param = nu)

test_that("`JoeExtArch2FParam`-class is correctly initialized", {
  parm <- JoeExtArch2FParam()
  expect_s4_class(parm, "JoeExtArch2FParam")

  parm@survival <- TRUE
  parm@copula <- copula::archmCopula(family = "Joe")
  setDimension(parm) <- d
  setLambda(parm) <- lambda
  setNu(parm) <- nu
  expect_true(validObject(parm))
  expect_equal(getDimension(parm), d)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getTau(parm), tau)

  expect_equal(parm, JoeExtArch2FParam(d, lambda, nu))
  expect_equal(parm, JoeExtArch2FParam(d, lambda, tau = tau))
})
