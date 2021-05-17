d <- 4L
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
  expect_error(validObject(parm), NA)
  expect_equal(getDimension(parm), d)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getTau(parm), tau)

  expect_equal(parm, JoeExtArch2FParam(d, lambda, nu))
  expect_equal(parm, JoeExtArch2FParam(d, lambda, tau = tau))
})

test_that("`JoeExtArch2FParam`-class setters can be used in arbitrary order", { # nolint
  parm <- JoeExtArch2FParam(d, lambda, nu)

  parm2 <- JoeExtArch2FParam()
  setDimension(parm2) <- d
  setLambda(parm2) <- lambda
  setTau(parm2) <- tau
  expect_equal(parm, parm2)

  parm2 <- JoeExtArch2FParam()
  setDimension(parm2) <- d
  setNu(parm2) <- nu
  setLambda(parm2) <- lambda
  expect_equal(parm, parm2)

  parm2 <- JoeExtArch2FParam()
  setLambda(parm2) <- lambda
  setDimension(parm2) <- d
  setNu(parm2) <- nu
  expect_equal(parm, parm2)

  parm2 <- JoeExtArch2FParam()
  setNu(parm2) <- nu
  setDimension(parm2) <- d
  setLambda(parm2) <- lambda
  expect_equal(parm, parm2)

  parm2 <- JoeExtArch2FParam()
  setLambda(parm2) <- lambda
  setNu(parm2) <- nu
  setDimension(parm2) <- d
  expect_equal(parm, parm2)

  parm2 <- JoeExtArch2FParam()
  setNu(parm2) <- nu
  setLambda(parm2) <- lambda
  setDimension(parm2) <- d
  expect_equal(parm, parm2)
})
