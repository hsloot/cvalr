partition <- list(1L:2L, 3L:6L, 7L:8L)
composition <- purrr::map_int(partition, length)
lambda <- 8e-2
nu <- c(0.3, 0.8)

d <- sum(composition)
rho <- (6 / pi) * asin(nu / 2)
tau <- (2 / pi) * asin(nu)

test_that("`H2ExtGaussian3FParam`-class is correctly initialized", {
  parm <- H2ExtGaussian3FParam()
  expect_s4_class(parm, "H2ExtGaussian3FParam")

  setComposition(parm) <- composition
  setLambda(parm) <- lambda
  setNu(parm) <- nu
  expect_error(validObject(parm), NA)
  expect_equal(getDimension(parm), d)
  expect_equal(getComposition(parm), composition)
  expect_equal(getPartition(parm), partition)
  expect_equal(getLambda(parm), lambda)
  expect_equal(getNu(parm), nu)
  expect_equal(getRho(parm), rho)
  expect_equal(getTau(parm), tau)

  expect_equal(parm, H2ExtGaussian3FParam(composition, lambda, nu))
  expect_equal(parm, H2ExtGaussian3FParam(composition, lambda, rho = rho))
  expect_equal(parm, H2ExtGaussian3FParam(composition, lambda, tau = tau))
})

test_that("`H2ExtGaussian3FParam`-class setters can be used in arbitrary order", { # nolint
  parm <- H2ExtGaussian3FParam(composition, lambda, nu, fraction)

  parm2 <- H2ExtGaussian3FParam()
  setComposition(parm2) <- composition
  setLambda(parm2) <- lambda
  setRho(parm2) <- rho
  expect_equal(parm, parm2)

  parm2 <- H2ExtGaussian3FParam()
  setComposition(parm2) <- composition
  setLambda(parm2) <- lambda
  setTau(parm2) <- tau
  expect_equal(parm, parm2)

  parm2 <- H2ExtGaussian3FParam()
  setComposition(parm2) <- composition
  setRho(parm2) <- rho
  setLambda(parm2) <- lambda
  expect_equal(parm, parm2)

  parm2 <- H2ExtGaussian3FParam()
  setLambda(parm2) <- lambda
  setComposition(parm2) <- composition
  setRho(parm2) <- rho
  expect_equal(parm, parm2)

  parm2 <- H2ExtGaussian3FParam()
  setRho(parm2) <- rho
  setComposition(parm2) <- composition
  setLambda(parm2) <- lambda
  expect_equal(parm, parm2)

  parm2 <- H2ExtGaussian3FParam()
  setLambda(parm2) <- lambda
  setRho(parm2) <- rho
  setComposition(parm2) <- composition
  expect_equal(parm, parm2)

  parm2 <- H2ExtGaussian3FParam()
  setRho(parm2) <- rho
  setLambda(parm2) <- lambda
  setComposition(parm2) <- composition
  expect_equal(parm, parm2)
})

test_that("`expected_pcds_loss` works as expected for `H2ExtGaussian3FParam", {
  # HELPER START
  epcdslfn <- function(parm, times, recovery_rate) {
    qassert(times, "N+[0,)")
    qassert(recovery_rate, "N1[0,1]")

    (1 - recovery_rate) * pexp(times, rate = getLambda(parm))
  }
  # HELPER END

  parm <- H2ExtGaussian3FParam(composition = composition, lambda = lambda, nu = nu)
  times <- seq(25e-2, 5L, by = 25e-2)
  recovery_rate <- 0.4

  x <- expected_pcds_loss(parm, times, recovery_rate = recovery_rate)
  expect_numeric(x, any.missing = FALSE, lower = 0, upper = 1,
    len = length(times), sorted = TRUE)
  expect_equal(x, epcdslfn(parm, times, recovery_rate))
})
