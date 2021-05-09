d <- 5L
lambda <- 8e-2
alpha <- 4e-1

rho <- 3 * alpha / (4 - alpha)
tau <- alpha / (2 - alpha)

nu <- log2(2 - alpha)
bf <- ScaledBernsteinFunction(
  scale = lambda,
  original = AlphaStableBernsteinFunction(alpha = nu))
ex_qmatrix <- rmo::exQMatrix(bf, d)
ex_intensities <- rmo::exIntensities(bf, d = d)

test_that("`AlphaStableExtMO2FParam`-class is correctly initialized", {
  parm <- AlphaStableExtMO2FParam()
  expect_s4_class(parm, "AlphaStableExtMO2FParam")

  setDimension(parm) <- d
  setBernsteinFunction(parm) <- bf
  setLambda(parm) <- lambda
  setNu(parm) <- nu
  expect_true(validObject(parm))
  expect_equal(getDimension(parm), d)
  expect_equal(getExQMatrix(parm), ex_qmatrix)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getBernsteinFunction(parm), bf)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)
  expect_equal(getAlpha(parm), alpha)

  expect_equal(parm, AlphaStableExtMO2FParam(d, lambda, nu))
  expect_equal(parm, AlphaStableExtMO2FParam(d, lambda, rho = rho))
  expect_equal(parm, AlphaStableExtMO2FParam(d, lambda, tau = tau))
  expect_equal(parm, AlphaStableExtMO2FParam(d, lambda, alpha = alpha))
  expect_equal(as(parm, "ExtMOParam"), ExtMOParam(d, bf))
  expect_equal(as(parm, "ExMOParam"), ExMOParam(ex_intensities))
  expect_equal(as(parm, "ExMarkovParam"), ExMarkovParam(ex_qmatrix))
})