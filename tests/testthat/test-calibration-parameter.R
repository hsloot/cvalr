times <- seq(0.25, 1, by = 0.25)
lambda <- 0.08
alpha <- 0.4
rho <- 3 * alpha / (4 - alpha)
tau <- alpha / (2 - alpha)

ex_intensities <- lambda * c(1 - alpha, alpha)
qmatrix <- matrix(
  c(
    - (2*ex_intensities[[1]]+ex_intensities[[2]]), 2*ex_intensities[[1]], ex_intensities[[2]],
    0, - (ex_intensities[[1]]+ex_intensities[[2]]), (ex_intensities[[1]]+ex_intensities[[2]]),
    0, 0, 0
  ),
  nrow = 3, ncol = 3, byrow = TRUE
)

test_that("Biv. ExMarkovParam is initialized correctly", {
  parm <- ExMarkovParam(qmatrix = qmatrix)
  expect_equal(getDimension(parm), 2L)
  expect_equal(getQMatrix(parm), qmatrix)

  expect_equal(
    sapply(times, function(t) expected_loss(parm, t, function(x) x)),
    pexp(times, rate = lambda))
})

test_that("Biv. ExMOParam is initialized correctly", {
  parm <- ExMOParam(ex_intensities = ex_intensities)
  expect_equal(getDimension(parm), 2L)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  expect_equal(
    sapply(times, function(t) expected_loss(parm, t, function(x) x)),
    pexp(times, rate = lambda))
})

test_that("Biv. CuadrasAugeExtMO2FParam is initialized correctly", {
  nu <- alpha
  parm <- CuadrasAugeExtMO2FParam(dim = 2, lambda = lambda, nu)
  expect_equal(getDimension(parm), 2L)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)
  expect_equal(getAlpha(parm), alpha)

  setRho(parm) <- rho
  expect_equal(getNu(parm), nu)

  setTau(parm) <- tau
  expect_equal(getNu(parm), nu)

  setAlpha(parm) <- alpha
  expect_equal(getNu(parm), nu)

  setNu(parm) <- nu
  expect_equal(getBernsteinFunction(parm),
    ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = 1 - nu),
        second = ConstantBernsteinFunction(constant = nu))))
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  rho <- 3 * alpha / (4 - alpha)
  parm <- CuadrasAugeExtMO2FParam(dim = 2, lambda = lambda, rho = rho)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  tau <- alpha / (2 - alpha)
  parm <- CuadrasAugeExtMO2FParam(dim = 2, lambda = lambda, tau = tau)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  parm <- CuadrasAugeExtMO2FParam(dim = 2, lambda = lambda, alpha = alpha)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  expect_equal(
    sapply(times, function(t) expected_loss(parm, t, function(x) x)),
    pexp(times, rate = lambda))
})

test_that("Biv. AlphaStableExtMO2FParam is initialized correctly", {
  nu <- log2(2 - alpha)
  parm <- AlphaStableExtMO2FParam(dim = 2, lambda = lambda, nu)
  expect_equal(getDimension(parm), 2L)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)
  expect_equal(getAlpha(parm), alpha)

  setRho(parm) <- rho
  expect_equal(getNu(parm), nu)

  setTau(parm) <- tau
  expect_equal(getNu(parm), nu)

  setAlpha(parm) <- alpha
  expect_equal(getNu(parm), nu)

  setNu(parm) <- nu
  expect_equal(getBernsteinFunction(parm),
               ScaledBernsteinFunction(
                 scale = lambda,
                 original = AlphaStableBernsteinFunction(alpha = nu)))
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  rho <- 3 * alpha / (4 - alpha)
  parm <- AlphaStableExtMO2FParam(dim = 2, lambda = lambda, rho = rho)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  tau <- alpha / (2 - alpha)
  parm <- AlphaStableExtMO2FParam(dim = 2, lambda = lambda, tau = tau)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  parm <- AlphaStableExtMO2FParam(dim = 2, lambda = lambda, alpha = alpha)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  expect_equal(
    sapply(times, function(t) expected_loss(parm, t, function(x) x)),
    pexp(times, rate = lambda))
})

test_that("Biv. PoissonExtMO2FParam is initialized correctly", {
  nu <- -log(1 - sqrt(alpha))
  parm <- PoissonExtMO2FParam(dim = 2, lambda = lambda, nu)
  expect_equal(getDimension(parm), 2L)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)
  expect_equal(getAlpha(parm), alpha)

  setRho(parm) <- rho
  expect_equal(getNu(parm), nu)

  setTau(parm) <- tau
  expect_equal(getNu(parm), nu)

  setAlpha(parm) <- alpha
  expect_equal(getNu(parm), nu)

  setNu(parm) <- nu
  expect_equal(getBernsteinFunction(parm),
               ScaledBernsteinFunction(
                 scale = lambda,
                 original = SumOfBernsteinFunctions(
                   first = LinearBernsteinFunction(scale = exp(-nu)),
                   second = PoissonBernsteinFunction(lambda = 1, eta = nu))))
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  rho <- 3 * alpha / (4 - alpha)
  parm <- PoissonExtMO2FParam(dim = 2, lambda = lambda, rho = rho)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  tau <- alpha / (2 - alpha)
  parm <- PoissonExtMO2FParam(dim = 2, lambda = lambda, tau = tau)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  parm <- PoissonExtMO2FParam(dim = 2, lambda = lambda, alpha = alpha)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  expect_equal(
    sapply(times, function(t) expected_loss(parm, t, function(x) x)),
    pexp(times, rate = lambda))
})

test_that("Biv. ExponentialExtMO2FParam is initialized correctly", {
  nu <- 0.5 * (-3 + sqrt(1 + 8 / alpha))
  parm <- ExponentialExtMO2FParam(dim = 2, lambda = lambda, nu)
  expect_equal(getDimension(parm), 2L)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)
  expect_equal(getAlpha(parm), alpha)

  setRho(parm) <- rho
  expect_equal(getNu(parm), nu)

  setTau(parm) <- tau
  expect_equal(getNu(parm), nu)

  setAlpha(parm) <- alpha
  expect_equal(getNu(parm), nu)

  setNu(parm) <- nu
  expect_equal(getBernsteinFunction(parm),
               ScaledBernsteinFunction(
                 scale = lambda,
                 original = SumOfBernsteinFunctions(
                   first = LinearBernsteinFunction(scale = 1 - 1 / (1 + nu)),
                   second = ExponentialBernsteinFunction(lambda = nu))))
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  rho <- 3 * alpha / (4 - alpha)
  parm <- ExponentialExtMO2FParam(dim = 2, lambda = lambda, rho = rho)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  tau <- alpha / (2 - alpha)
  parm <- ExponentialExtMO2FParam(dim = 2, lambda = lambda, tau = tau)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  parm <- ExponentialExtMO2FParam(dim = 2, lambda = lambda, alpha = alpha)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  expect_equal(
    sapply(times, function(t) expected_loss(parm, t, function(x) x)),
    pexp(times, rate = lambda))
})
