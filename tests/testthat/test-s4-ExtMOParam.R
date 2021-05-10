d <- 5L
lambda <- 8e-2
alpha <- 4e-1

nu <- log2(2 - alpha)
bf <- ScaledBernsteinFunction(
  scale = lambda, original = AlphaStableBernsteinFunction(alpha = nu))
ex_qmatrix <- rmo::exQMatrix(bf, d)
ex_intensities <- rmo::exIntensities(bf, d = d)

test_that("`ExtMOParam`-class is correctly initialized", {
  parm <- ExtMOParam()
  expect_s4_class(parm, "ExtMOParam")

  setDimension(parm) <- d
  setBernsteinFunction(parm) <- bf
  expect_true(validObject(parm))
  expect_equal(getDimension(parm), d)
  expect_equal(getExQMatrix(parm), ex_qmatrix)
  expect_equal(getExIntensities(parm), ex_intensities)
  expect_equal(getBernsteinFunction(parm), bf)

  expect_equal(parm, ExtMOParam(d, bf))
  expect_equal(as(parm, "ExMOParam"), ExMOParam(ex_intensities))
  expect_equal(as(parm, "ExMarkovParam"), ExMarkovParam(ex_qmatrix))
})
