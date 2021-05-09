lambda <- 0.08
alpha <- 0.4

rho <- 3 * alpha / (4 - alpha)
tau <- alpha / (2 - alpha)

times <- seq(0.25, 1, by = 0.25)

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
