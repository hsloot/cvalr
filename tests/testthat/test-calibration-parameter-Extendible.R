lambda <- 0.08
alpha <- 0.4

rho <- 3 * alpha / (4 - alpha)
tau <- alpha / (2 - alpha)
ex_intensities <- lambda * c(2 * (1 - alpha), alpha)
ex_qmatrix <- matrix(
  c(
    - (ex_intensities[[1]]+ex_intensities[[2]]), ex_intensities[[1]], ex_intensities[[2]],
    0, - (0.5 * ex_intensities[[1]] + ex_intensities[[2]]),
      (0.5 * ex_intensities[[1]] + ex_intensities[[2]]),
    0, 0, 0
  ),
  nrow = 3, ncol = 3, byrow = TRUE
)

times <- seq(0.25, 1, by = 0.25)

test_that("Biv. CuadrasAugeExtMO2FParam is initialized correctly", {
  nu <- alpha
  parm <- CuadrasAugeExtMO2FParam(dim = 2L, lambda = 5e-2, rho = 5e-1)
  setLambda(parm) <- lambda
  setNu(parm) <- nu
  expect_equal(getDimension(parm), 2L)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)
  expect_equal(getAlpha(parm), alpha)

  setAlpha(parm) <- alpha
  expect_equal(getNu(parm), nu)

  setRho(parm) <- rho
  expect_equal(getNu(parm), nu)

  setTau(parm) <- tau
  expect_equal(getNu(parm), nu)

  expect_equal(getBernsteinFunction(parm),
    ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = 1 - nu),
        second = ConstantBernsteinFunction(constant = nu))))
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)

  parm <- CuadrasAugeExtMO2FParam(dim = 2, lambda = lambda, rho = rho)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)

  parm <- CuadrasAugeExtMO2FParam(dim = 2, lambda = lambda, tau = tau)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)

  parm <- CuadrasAugeExtMO2FParam(dim = 2, lambda = lambda, alpha = alpha)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)
})

test_that("Biv. AlphaStableExtMO2FParam is initialized correctly", {
  nu <- log2(2 - alpha)
  parm <- AlphaStableExtMO2FParam(dim = 2L, lambda = 5e-2, rho = 5e-1)
  setLambda(parm) <- lambda
  setNu(parm) <- nu
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

  expect_equal(getBernsteinFunction(parm),
               ScaledBernsteinFunction(
                 scale = lambda,
                 original = AlphaStableBernsteinFunction(alpha = nu)))
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)

  parm <- AlphaStableExtMO2FParam(dim = 2L, lambda = lambda, rho = rho)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)

  parm <- AlphaStableExtMO2FParam(dim = 2L, lambda = lambda, tau = tau)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)

  parm <- AlphaStableExtMO2FParam(dim = 2L, lambda = lambda, alpha = alpha)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)
})

test_that("Biv. PoissonExtMO2FParam is initialized correctly", {
  nu <- -log(1 - sqrt(alpha))
  parm <- PoissonExtMO2FParam(dim = 2L, lambda = 5e-2, rho = 5e-1)
  setLambda(parm) <- lambda
  setNu(parm) <- nu
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

  expect_equal(getBernsteinFunction(parm),
               ScaledBernsteinFunction(
                 scale = lambda,
                 original = SumOfBernsteinFunctions(
                   first = LinearBernsteinFunction(scale = exp(-nu)),
                   second = PoissonBernsteinFunction(lambda = 1, eta = nu))))
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)

  parm <- PoissonExtMO2FParam(dim = 2L, lambda = lambda, rho = rho)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)

  parm <- PoissonExtMO2FParam(dim = 2L, lambda = lambda, tau = tau)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)

  parm <- PoissonExtMO2FParam(dim = 2L, lambda = lambda, alpha = alpha)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)
})

test_that("Biv. ExponentialExtMO2FParam is initialized correctly", {
  nu <- 0.5 * (-3 + sqrt(1 + 8 / alpha))
  parm <- ExponentialExtMO2FParam(dim = 2L, lambda = 5e-2, rho = 5e-1)
  setLambda(parm) <- lambda
  setNu(parm) <- nu
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

  expect_equal(getBernsteinFunction(parm),
               ScaledBernsteinFunction(
                 scale = lambda,
                 original = SumOfBernsteinFunctions(
                   first = LinearBernsteinFunction(scale = 1 - 1 / (1 + nu)),
                   second = ExponentialBernsteinFunction(lambda = nu))))
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)

  parm <- ExponentialExtMO2FParam(dim = 2L, lambda = lambda, rho = rho)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)

  parm <- ExponentialExtMO2FParam(dim = 2L, lambda = lambda, tau = tau)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)

  parm <- ExponentialExtMO2FParam(dim = 2L, lambda = lambda, alpha = alpha)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getExQMatrix(parm), ex_qmatrix)
})

test_that("ExtGaussian2FParam is initialized correctly", {
  nu <- 2 * sin(pi / 6 * rho)
  tau <- 2 / pi * asin(nu)
  parm <- ExtGaussian2FParam(dim = 125L, lambda = lambda, nu)
  parm <- ExtGaussian2FParam(dim = 125L, lambda = lambda, rho = rho)
  parm <- ExtGaussian2FParam(dim = 125L, lambda = lambda, tau = tau)
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)

  setTau(parm) <- tau
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)

  setRho(parm) <- rho
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)
})


test_that("ClaytonExtArch2FParam is initialized correctly", {
  parm <- ClaytonExtArch2FParam(dim = 125L, lambda = lambda, rho = rho)
  tau <- getTau(parm)
  nu <- getNu(parm)
  parm <- ClaytonExtArch2FParam(dim = 125L, lambda = lambda, nu)
  parm <- ClaytonExtArch2FParam(dim = 125L, lambda = lambda, tau = tau)
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)

  setTau(parm) <- tau
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)

  setRho(parm) <- rho
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)
})

test_that("FrankExtArch2FParam is initialized correctly", {
  parm <- FrankExtArch2FParam(dim = 125L, lambda = lambda, rho = rho)
  tau <- getTau(parm)
  nu <- getNu(parm)
  parm <- FrankExtArch2FParam(dim = 125L, lambda = lambda, nu)
  parm <- FrankExtArch2FParam(dim = 125L, lambda = lambda, tau = tau)
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)

  setTau(parm) <- tau
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)

  setRho(parm) <- rho
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)
})

test_that("GumbelExtArch2FParam is initialized correctly", {
  parm <- GumbelExtArch2FParam(dim = 125L, lambda = lambda, rho = rho)
  tau <- getTau(parm)
  nu <- getNu(parm)
  parm <- GumbelExtArch2FParam(dim = 125L, lambda = lambda, nu)
  parm <- GumbelExtArch2FParam(dim = 125L, lambda = lambda, tau = tau)
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)

  setTau(parm) <- tau
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)

  setRho(parm) <- rho
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)
})

test_that("AmhExtArch2FParam is initialized correctly", {
  parm <- AmhExtArch2FParam(dim = 125L, lambda = lambda, rho = rho)
  tau <- getTau(parm)
  nu <- getNu(parm)
  parm <- AmhExtArch2FParam(dim = 125L, lambda = lambda, nu)
  parm <- AmhExtArch2FParam(dim = 125L, lambda = lambda, tau = tau)
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)

  setTau(parm) <- tau
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)

  setRho(parm) <- rho
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)
})

test_that("JoeExtArch2FParam is initialized correctly", {
  parm <- JoeExtArch2FParam(dim = 125L, lambda = lambda, tau = tau)
  nu <- getNu(parm)
  parm <- JoeExtArch2FParam(dim = 125L, lambda = lambda, nu)
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)

  setTau(parm) <- tau
  expect_equal(getDimension(parm), 125L)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)
})
