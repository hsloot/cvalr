lambda <- 0.08
alpha <- 0.4

rho <- 3 * alpha / (4 - alpha)
tau <- alpha / (2 - alpha)
ex_intensities <- lambda * c(2 * (1 - alpha), alpha)
qmatrix <- matrix(
  c(
    - (ex_intensities[[1]]+ex_intensities[[2]]), ex_intensities[[1]], ex_intensities[[2]],
    0, - (0.5 * ex_intensities[[1]] + ex_intensities[[2]]), (0.5 * ex_intensities[[1]] + ex_intensities[[2]]),
    0, 0, 0
  ),
  nrow = 3, ncol = 3, byrow = TRUE
)

times <- seq(0.25, 1, by = 0.25)

test_that("Biv. ExMarkovParam is initialized correctly", {
  parm <- ExMarkovParam(qmatrix = qmatrix)
  expect_equal(getDimension(parm), 2L)
  expect_equal(getQMatrix(parm), qmatrix)

  expect_equal(
    sapply(times,
      function(t) {
        expected_pcds_loss(parm, t, recovery_rate = 0)
        }),
    pexp(times, rate = lambda))

  expect_equal(
    expected_pcds_equation(parm,
      times = times, discount_factors = rep_len(1, length(times)),
      recovery_rate = 0.4, coupon = 0.08, upfront = 0.01),
    portfolio_cds_equation(
      expected_losses = expected_pcds_loss(parm, times = times, recovery_rate=0.4),
      times = times, discount_factors = rep_len(1, length(times)),
      recovery_rate = 0.4, coupon = 0.08, upfront = 0.01)
  )
  expect_equal(
    expected_cdo_equation(parm,
      times = times, discount_factors = rep_len(1, length(times)),
      recovery_rate = 0.4, lower = 0.1, upper = 0.2,
      coupon = 0.08, upfront = 0.01),
    cdo_equation(
      expected_losses = expected_cdo_loss(parm, times = times,
        recovery_rate = 0.4, lower = 0.1, upper = 0.2),
      times=times, discount_factors = rep_len(1, length(times)),
      lower = 0.1, upper = 0.2, coupon = 0.08, upfront = 0.01)
  )
})

test_that("Biv. ExMOParam is initialized correctly", {
  parm <- ExMOParam(ex_intensities = ex_intensities)
  expect_equal(getDimension(parm), 2L)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getQMatrix(parm), qmatrix)

  expect_equal(
    expected_pcds_loss(parm, times, recovery_rate = 0),
    pexp(times, rate = lambda))
})

test_that("Biv. CuadrasAugeExtMO2FParam is initialized correctly", {
  nu <- alpha
  parm <- CuadrasAugeExtMO2FParam(dim = 2)
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
    expected_pcds_loss(parm, times, recovery_rate = 0),
    pexp(times, rate = lambda))
  expect_equal(
    expected_pcds_loss(parm, times, recovery_rate = 0, method = "CalibrationParam"),
    pexp(times, rate = lambda))
})

test_that("Biv. AlphaStableExtMO2FParam is initialized correctly", {
  nu <- log2(2 - alpha)
  parm <- AlphaStableExtMO2FParam(dim = 2)
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
    expected_pcds_loss(parm, times, recovery_rate = 0),
    pexp(times, rate = lambda))
  expect_equal(
    expected_pcds_loss(parm, times, recovery_rate = 0, method = "CalibrationParam"),
    pexp(times, rate = lambda))
})

test_that("Biv. PoissonExtMO2FParam is initialized correctly", {
  nu <- -log(1 - sqrt(alpha))
  parm <- PoissonExtMO2FParam(dim = 2)
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
    expected_pcds_loss(parm, times, recovery_rate = 0),
    pexp(times, rate = lambda))
  expect_equal(
    expected_pcds_loss(parm, times, recovery_rate = 0, method = "CalibrationParam"),
    pexp(times, rate = lambda))
})

test_that("Biv. ExponentialExtMO2FParam is initialized correctly", {
  nu <- 0.5 * (-3 + sqrt(1 + 8 / alpha))
  parm <- ExponentialExtMO2FParam(dim = 2)
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
    expected_pcds_loss(parm, times, recovery_rate = 0),
    pexp(times, rate = lambda))
  expect_equal(
    expected_pcds_loss(parm, times, recovery_rate = 0, method = "CalibrationParam"),
    pexp(times, rate = lambda))
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

  expect_equal(
    expected_pcds_loss(parm, times, recovery_rate = 0),
    pexp(times, rate = lambda),
    tolerance = 1e-2)
  expect_equal(
    expected_pcds_loss(parm, times, recovery_rate = 0, method = "CalibrationParam"),
    pexp(times, rate = lambda),
    tolerance = 1e-2)

  recovery_rate <- 0.4
  lower <- 0.1
  upper <- 0.2
  expect_equal(
    expected_value(parm, times, function(x) {
      pmin(pmax((1 - recovery_rate) * x - lower, 0), upper - lower)
    }),
    expected_cdo_loss(parm, times, recovery_rate, lower, upper, method = "default"),
    tolerance = 1e-2
  )
  expect_equal(
    expected_value(parm, times, function(x) {
      pmin(pmax((1 - recovery_rate) * x - lower, 0), upper - lower)
    }),
    expected_cdo_loss(parm, times, recovery_rate, lower, upper, method = "CalibrationParam"),
    tolerance = 1e-2
  )
})


test_that("FrankExtArch2FParam is initialized correctly", {
  parm <- FrankExtArch2FParam(dim = 125L, lambda = lambda, rho = rho)
  tau <- getTau(parm)
  nu <- getNu(parm)
  parm <- FrankExtArch2FParam(dim = 125L, lambda = lambda, nu)
  parm <- FrankExtArch2FParam(dim = 125L, lambda = lambda, tau = tau)
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

  expect_equal(
    expected_pcds_loss(parm, times, recovery_rate = 0),
    pexp(times, rate = lambda),
    tolerance = 1e-2)
  expect_equal(
    expected_pcds_loss(parm, times, recovery_rate = 0, pd_args = list(n_sim = 1e4, seed = 1623)),
    pexp(times, rate = lambda),
    tolerance = 1e-2)

  recovery_rate <- 0.4
  lower <- 0.1
  upper <- 0.2
  expect_equal(
    expected_value(parm, times, function(x) {
      pmin(pmax((1 - recovery_rate) * x - lower, 0), upper - lower)
    }, pd_args = list(n_sim = 1e4, seed = 1623)),
    expected_cdo_loss(parm, times, recovery_rate, lower, upper, pd_args = list(n_sim = 1e4, seed = 1623)),
    tolerance = 1e-2
  )
  expect_equal(
    expected_value(parm, times, function(x) {
      pmin(pmax((1 - recovery_rate) * x - lower, 0), upper - lower)
    }, pd_args = list(n_sim = 1e4, seed = 1623)),
    expected_cdo_loss(parm, times, recovery_rate, lower, upper, pd_args = list(n_sim = 1e4, seed = 1623)),
    tolerance = 1e-2
  )
})
