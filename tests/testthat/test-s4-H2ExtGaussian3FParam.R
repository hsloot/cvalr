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

test_that("`expected_pcds_equation` works as expected for `H2ExtGaussian3FParam", {
  times <- seq(0.25, 5, by = 0.25)
  discount_factors <- rep(1, length(times))
  recovery_rate <- 4e-1
  coupon <- 1e-1
  upfront <- -1e-2
  parm <- H2ExtGaussian3FParam(composition = composition, lambda = lambda, nu = nu)

  # using default
  x <- expected_pcds_equation(parm, times, discount_factors, recovery_rate, coupon, upfront)
  expect_numeric(x, finite = TRUE, any.missing = FALSE, len = 1L)
  y <- test__expected_pcds_equation__default(
    parm, times, discount_factors, recovery_rate, coupon, upfront)
  expect_equal(x, y)
})
