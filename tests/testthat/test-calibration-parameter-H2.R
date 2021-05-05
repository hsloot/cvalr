rep_last <- function(x, n = 1L) {
  c(x[-length(x)], rep_len(x[[length(x)]], length.out = n))
}

fraction <- 0.7
partition <- list(1:2, 3:6, 7:8)
lambda <- 8e-2
alpha0 <- c(0.3, 0.8)

dim <- length(unlist(partition))

alpha <- cumsum(c(fraction, 1 - fraction) * alpha0)
rho <- 3 * alpha / (4 - alpha)
tau <- alpha / (2 - alpha)

alpha0 <- rep_last(c(alpha[[1]], diff(alpha)) / c(fraction, 1 - fraction),
                    n = length(partition))
rho0 <- purrr::map_dbl(alpha0, ~{
  3 * .x / (4 - .x)
  })
tau0 <- purrr::map_dbl(alpha0, ~{
  .x / (2 - .x)
  })

ex_qmatrix_models <- purrr::map2(c(length(unlist(partition)), purrr::map_int(partition, length)),
                              alpha0, ~{
    as(new("AlphaStableExtMO2FParam", dim = .x, lambda = lambda, alpha = .y), "ExMarkovParam")
  })

test_that("H2ExMarkovParam is initialized correctly", {
  parm <- H2ExMarkovParam(models = ex_qmatrix_models, fraction = fraction)
  expect_equal(getDimension(parm), dim)
  expect_equal(getFraction(parm), fraction)
  expect_equal(getPartition(parm), partition)
  purrr::walk2(getModels(parm), ex_qmatrix_models, ~{
      expect_equal(getExQMatrix(.x), getExQMatrix(.y))
    })
})

ex_intensities_models <- purrr::map2(
  c(length(unlist(partition)), purrr::map_int(partition, length)),
  alpha0, ~{
    as(new("AlphaStableExtMO2FParam", dim = .x, lambda = lambda, alpha = .y), "ExMOParam")
  })

test_that("H2ExMOParam is initialized correctly", {
  parm <- H2ExMOParam(models = ex_intensities_models, fraction = fraction)
  expect_equal(getDimension(parm), dim)
  expect_equal(getFraction(parm), fraction)
  expect_equal(getPartition(parm), partition)
  purrr::walk2(getModels(parm), ex_intensities_models, ~{
      expect_equal(getExIntensities(.x), getExIntensities(.y))
    })
})

ext_bf_models <- purrr::map2(c(length(unlist(partition)), purrr::map_int(partition, length)),
                              alpha0, ~{
    as(new("AlphaStableExtMO2FParam", dim = .x, lambda = lambda, alpha = .y), "ExtMOParam")
  })

test_that("H2ExtMOParam is initialized correctly", {
  parm <- H2ExtMOParam(models = ext_bf_models, fraction = fraction)
  expect_equal(getDimension(parm), dim)
  expect_equal(getFraction(parm), fraction)
  expect_equal(getPartition(parm), partition)
  purrr::walk2(getModels(parm), ext_bf_models, ~{
      expect_equal(getBernsteinFunction(.x), getBernsteinFunction(.y))
    })
})

test_that("CuadrasAugeH2ExtMO3FParam is initialized correctly", {
  nu <- c(alpha[[1]], diff(alpha)) / c(fraction, 1 - fraction)
  parm <- CuadrasAugeH2ExtMO3FParam(partition = partition, fraction = fraction,
                                    lambda = lambda, nu = nu)
  expect_equal(getDimension(parm), dim)
  expect_equal(getFraction(parm), fraction)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)
  expect_equal(getAlpha(parm), alpha)

  expect_equal(nu, getNu(
    CuadrasAugeH2ExtMO3FParam(partition = partition, fraction = fraction,
                              lambda = lambda, rho = rho)))
  expect_equal(nu, getNu(
    CuadrasAugeH2ExtMO3FParam(partition = partition, fraction = fraction,
                              lambda = lambda, tau = tau)))
  expect_equal(nu, getNu(
    CuadrasAugeH2ExtMO3FParam(partition = partition, fraction = fraction,
                              lambda = lambda, alpha = alpha)))
})

test_that("AlphaStableH2ExtMO3FParam is initialized correctly", {
  nu <- log2(2 - c(alpha[[1]], diff(alpha)) / c(fraction, 1 - fraction))
  parm <- AlphaStableH2ExtMO3FParam(partition = partition, fraction = fraction,
                                    lambda = lambda, nu = nu)
  expect_equal(getDimension(parm), dim)
  expect_equal(getFraction(parm), fraction)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)
  expect_equal(getAlpha(parm), alpha)

  expect_equal(nu, getNu(
    AlphaStableH2ExtMO3FParam(partition = partition, fraction = fraction,
                              lambda = lambda, rho = rho)))
  expect_equal(nu, getNu(
    AlphaStableH2ExtMO3FParam(partition = partition, fraction = fraction,
                              lambda = lambda, tau = tau)))
  expect_equal(nu, getNu(
    AlphaStableH2ExtMO3FParam(partition = partition, fraction = fraction,
                              lambda = lambda, alpha = alpha)))
})

test_that("PoissonH2ExtMO3FParam is initialized correctly", {
  nu <- -log(1 - sqrt(c(alpha[[1]], diff(alpha)) / c(fraction, 1 - fraction)))
  parm <- PoissonH2ExtMO3FParam(partition = partition, fraction = fraction,
                                    lambda = lambda, nu = nu)
  expect_equal(getDimension(parm), dim)
  expect_equal(getFraction(parm), fraction)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)
  expect_equal(getAlpha(parm), alpha)

  expect_equal(nu, getNu(
    PoissonH2ExtMO3FParam(partition = partition, fraction = fraction,
                              lambda = lambda, rho = rho)))
  expect_equal(nu, getNu(
    PoissonH2ExtMO3FParam(partition = partition, fraction = fraction,
                              lambda = lambda, tau = tau)))
  expect_equal(nu, getNu(
    PoissonH2ExtMO3FParam(partition = partition, fraction = fraction,
                              lambda = lambda, alpha = alpha)))
})

test_that("ExponentialH2ExtMO3FParam is initialized correctly", {
  nu <- 0.5 * (-3 + sqrt(1 + 8 / (c(alpha[[1]], diff(alpha)) / c(fraction, 1 - fraction))))
  parm <- ExponentialH2ExtMO3FParam(partition = partition, fraction = fraction,
                                    lambda = lambda, nu = nu)
  expect_equal(getDimension(parm), dim)
  expect_equal(getFraction(parm), fraction)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)
  expect_equal(getAlpha(parm), alpha)

  expect_equal(nu, getNu(
    ExponentialH2ExtMO3FParam(partition = partition, fraction = fraction,
                          lambda = lambda, rho = rho)))
  expect_equal(nu, getNu(
    ExponentialH2ExtMO3FParam(partition = partition, fraction = fraction,
                          lambda = lambda, tau = tau)))
  expect_equal(nu, getNu(
    ExponentialH2ExtMO3FParam(partition = partition, fraction = fraction,
                          lambda = lambda, alpha = alpha)))
})

test_that("H2ExtGaussian3FParam is initialized correctly", {
  nu <- 2 * sin(pi / 6 * rho)
  tau <- 2 / pi * asin(nu)
  parm <- H2ExtGaussian3FParam(partition = partition,
                               lambda = lambda, nu = nu)
  expect_equal(getDimension(parm), dim)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)

  expect_equal(nu, getNu(
    H2ExtGaussian3FParam(partition = partition, lambda = lambda, rho = rho)))
  expect_equal(nu, getNu(
    H2ExtGaussian3FParam(partition = partition, lambda = lambda, tau = tau)))
})

test_that("ClaytonH2ExtArch3FParam is initialized correctly", {
  parm <- ClaytonH2ExtArch3FParam(partition = partition,
                                lambda = lambda, rho = rho)
  nu <- getNu(parm)
  tau <- getTau(parm)
  parm <- ClaytonH2ExtArch3FParam(partition = partition,
                                lambda = lambda, nu = nu)
  expect_equal(getDimension(parm), dim)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)

  expect_equal(nu, getNu(
    ClaytonH2ExtArch3FParam(partition = partition, lambda = lambda, rho = rho)))
  expect_equal(nu, getNu(
    ClaytonH2ExtArch3FParam(partition = partition, lambda = lambda, tau = tau)))
})

test_that("FrankH2ExtArch3FParam is initialized correctly", {
  parm <- FrankH2ExtArch3FParam(partition = partition,
                                lambda = lambda, rho = rho)
  nu <- getNu(parm)
  tau <- getTau(parm)
  parm <- FrankH2ExtArch3FParam(partition = partition,
                                lambda = lambda, nu = nu)
  expect_equal(getDimension(parm), dim)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)

  expect_equal(nu, getNu(
    FrankH2ExtArch3FParam(partition = partition, lambda = lambda, rho = rho)))
  expect_equal(nu, getNu(
    FrankH2ExtArch3FParam(partition = partition, lambda = lambda, tau = tau)))
})

test_that("GumbelH2ExtArch3FParam is initialized correctly", {
  parm <- GumbelH2ExtArch3FParam(partition = partition,
                                lambda = lambda, rho = rho)
  nu <- getNu(parm)
  tau <- getTau(parm)
  parm <- GumbelH2ExtArch3FParam(partition = partition,
                                lambda = lambda, nu = nu)
  expect_equal(getDimension(parm), dim)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)

  expect_equal(nu, getNu(
    GumbelH2ExtArch3FParam(partition = partition, lambda = lambda, rho = rho)))
  expect_equal(nu, getNu(
    GumbelH2ExtArch3FParam(partition = partition, lambda = lambda, tau = tau)))
})

test_that("AmhH2ExtArch3FParam is initialized correctly", {
  parm <- AmhH2ExtArch3FParam(partition = partition,
                                lambda = lambda, rho = rho)
  nu <- getNu(parm)
  tau <- getTau(parm)
  parm <- AmhH2ExtArch3FParam(partition = partition,
                                lambda = lambda, nu = nu)
  expect_equal(getDimension(parm), dim)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getRho(parm), rho, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)

  expect_equal(nu, getNu(
    AmhH2ExtArch3FParam(partition = partition, lambda = lambda, rho = rho)))
  expect_equal(nu, getNu(
    AmhH2ExtArch3FParam(partition = partition, lambda = lambda, tau = tau)))
})

test_that("JoeH2ExtArch3FParam is initialized correctly", {
  parm <- JoeH2ExtArch3FParam(partition = partition,
                                lambda = lambda, tau = tau)
  nu <- getNu(parm)
  tau <- getTau(parm)
  parm <- JoeH2ExtArch3FParam(partition = partition,
                                lambda = lambda, nu = nu)
  expect_equal(getDimension(parm), dim)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda, tolerance = 1e-2)
  expect_equal(getNu(parm), nu, tolerance = 1e-2)
  expect_equal(getTau(parm), tau, tolerance = 1e-2)

  expect_equal(nu, getNu(
    JoeH2ExtArch3FParam(partition = partition, lambda = lambda, tau = tau)))
})
