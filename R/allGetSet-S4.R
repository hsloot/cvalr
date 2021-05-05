#' @include allClass-S4.R allGeneric-S4.R checkmate.R
NULL

setMethod("getDimension", "CalibrationParam",
  function(object) {
    object@dim
  })

#' @importFrom checkmate qassert
setReplaceMethod("setDimension", "CalibrationParam",
  function(object, value) {
    qassert(value, "X1(0,)")
    object@dim <- as.integer(value)

    invisible(object)
  })
#' @importFrom checkmate qassert
setReplaceMethod("setDimension", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "X1(0,)")
    object@dim <- as.integer(value)
    object@copula@dimension <- as.integer(value)

    invisible(object)
  })

setMethod("getExQMatrix", "ExMarkovParam",
  function(object) {
    object@ex_qmatrix
  })

setReplaceMethod("setExQMatrix", "ExMarkovParam",
  function(object, value) {
    assert_exqmatrix(value, min.rows = 1L, min.cols = 1L)

    dim <- nrow(value)-1
    setDimension(object) <- dim
    object@ex_qmatrix <- value

    invisible(object)
  })


setMethod("getExIntensities", "ExMOParam",
  function(object) {
    object@ex_intensities
  })

#' @importFrom checkmate qassert
setReplaceMethod("setExIntensities", "ExMOParam",
  function(object, value) {
    qassert(value, "N+[0,)")
    qassert(max(value), "N1(0,)")
    setDimension(object) <- length(value)
    object@ex_intensities <- value
    setExQMatrix(object) <- rmo:::exi2exqm(value)

    invisible(object)
  })


setMethod("getBernsteinFunction", "ExtMOParam",
  function(object) {
    object@bf
  })

#' @importFrom rmo exIntensities
#' @importFrom checkmate assert_class
setReplaceMethod("setBernsteinFunction", "ExtMOParam",
  function(object, value) {
    assert_class(value, "BernsteinFunction")
    object@bf <- value
    setExIntensities(object) <- exIntensities(object@bf, object@dim)

    invisible(object)
  })
#' @importFrom rmo ScaledBernsteinFunction valueOf
#' @importFrom checkmate assert check_class
setReplaceMethod("setBernsteinFunction", "ExtMO2FParam",
  function(object, value) {
    assert(combine = "and",
      check_class(value, "ScaledBernsteinFunction"),
      check_equal(1, valueOf(value@original, 1, 0L)))
    object@lambda <- value@scale
    object@nu <- invAlpha(object, 2 - valueOf(value@original, 2, 0L))

    callNextMethod(object, value)
  })


setMethod("getLambda", "ExtMO2FParam",
  function(object) {
    object@lambda
  })
setMethod("getLambda", "ExtGaussian2FParam",
  function(object) {
    object@lambda
  })
setMethod("getLambda", "ExtArch2FParam",
  function(object) {
    object@lambda
  })
setMethod("getLambda", "H2ExtMO3FParam",
  function(object) {
    object@lambda
  })
setMethod("getLambda", "H2ExtGaussian3FParam",
  function(object) {
    object@lambda
  })
setMethod("getLambda", "H2ExtArch3FParam",
  function(object) {
    object@lambda
  })
#' @importFrom checkmate qassert
setReplaceMethod("setLambda", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    setBernsteinFunction(object) <- constructBernsteinFunction(
      object, value, object@nu)

    invisible(object)
  })
#' @importFrom checkmate qassert
setReplaceMethod("setLambda", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    object@lambda <- value

    invisible(object)
  })
#' @importFrom checkmate qassert
setReplaceMethod("setLambda", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    object@lambda <- value

    invisible(object)
  })


setMethod("getNu", "ExtMO2FParam",
  function(object) {
    object@nu
  })
setMethod("getNu", "ExtGaussian2FParam",
  function(object) {
    object@nu
  })
setMethod("getNu", "ExtArch2FParam",
  function(object) {
    object@nu
  })
setMethod("getNu", "H2ExtMO3FParam",
  function(object) {
    object@nu
  })
setMethod("getNu", "H2ExtGaussian3FParam",
  function(object) {
    object@nu
  })
setMethod("getNu", "H2ExtArch3FParam",
  function(object) {
    object@nu
  })

#' @importFrom checkmate qassert
setReplaceMethod("setNu", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1")
    setBernsteinFunction(object) <- constructBernsteinFunction(
      object, object@lambda, value)

    invisible(object)
  })
#' @importFrom checkmate qassert
setReplaceMethod("setNu", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    object@nu <- value

    invisible(object)
  })
#' @importFrom copula setTheta
#' @importFrom checkmate qassert
setReplaceMethod("setNu", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1")
    object@nu <- value
    object@copula <- setTheta(object@copula, value)

    invisible(object)
  })
#' @importFrom copula onacopulaL
#' @importFrom purrr map
#' @importFrom checkmate qassert
setReplaceMethod("setNu", "H2ExtArch3FParam",
  function(object, value) {
    qassert(value, "N2")
    object@nu <- value
    object@copula <- onacopulaL(
      family = object@family,
      nacList = list(value[[1]], map(object@partition, ~list(list(value[[2]], .)))))

    invisible(object)
  })


setMethod("getRho", "ExtMO2FParam",
  function(object) {
    alpha <- getAlpha(object)

    3 * alpha / (4 - alpha)
  })
setMethod("getRho", "ExtGaussian2FParam",
  function(object) {
    (6 / pi) * asin(getNu(object) / 2)
  })
#' @importFrom copula rho frankCopula
setMethod("getRho", "ExtArch2FParam",
  function(object) {
    copula::rho(object@copula)
  })
setMethod("getRho", "H2ExtMO3FParam",
  function(object) {
    alpha <- getAlpha(object)

    3 * alpha / (4 - alpha)
  })
setMethod("getRho", "H2ExtGaussian3FParam",
  function(object) {
    (6 / pi) * asin(getNu(object) / 2)
  })
#' @importFrom copula rho
#' @importFrom purrr map_dbl
setMethod("getRho", "H2ExtArch3FParam",
  function(object) {
    rho <- function(x) {
      copula::rho(archmCopula(family = object@family, param = x@theta, dim = 2))
    }
    c(rho(object@copula@copula), rho(object@copula@childCops[[1]]@copula))
  })

#' @importFrom checkmate qassert
setReplaceMethod("setRho", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invRho(object, value)

    invisible(object)
  })
#' @importFrom checkmate qassert
setReplaceMethod("setRho", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invRho(object, value)

    invisible(object)
  })
#' @importFrom checkmate qassert
setReplaceMethod("setRho", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invRho(object, value)

    invisible(object)
  })


setMethod("getTau", "ExtMO2FParam",
  function(object) {
    alpha <- getAlpha(object)

    alpha / (2 - alpha)
  })
setMethod("getTau", "ExtGaussian2FParam",
  function(object) {
    (2 / pi) * asin(getNu(object))
  })
#' @importFrom copula tau frankCopula
setMethod("getTau", "ExtArch2FParam",
  function(object) {
    copula::tau(object@copula)
  })
setMethod("getTau", "H2ExtMO3FParam",
  function(object) {
    alpha <- getAlpha(object)

    alpha / (2 - alpha)
  })
setMethod("getTau", "H2ExtGaussian3FParam",
  function(object) {
    (2 / pi) * asin(getNu(object))
  })
#' @importFrom copula rho
#' @importFrom purrr map_dbl
setMethod("getTau", "H2ExtArch3FParam",
  function(object) {
    tau <- function(x) {
      copula::tau(archmCopula(family = object@family, param = x@theta, dim = 2))
    }
    c(tau(object@copula@copula), tau(object@copula@childCops[[1]]@copula))
  })

#' @importFrom checkmate qassert
setReplaceMethod("setTau", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invTau(object, value)

    invisible(object)
  })
#' @importFrom checkmate qassert
setReplaceMethod("setTau", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invTau(object, value)

    invisible(object)
  })
#' @importFrom checkmate qassert
setReplaceMethod("setTau", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invTau(object, value)

    invisible(object)
  })


#' @importFrom rmo valueOf
setMethod("getAlpha", "ExtMO2FParam",
  function(object) {
    2 - valueOf(object@bf, 2, 0L) / valueOf(object@bf, 1, 0L)
  })
#' @importFrom purrr map_dbl
setMethod("getAlpha", "H2ExtMO3FParam",
  function(object) {
    fraction <- getFraction(object)
    alpha0 <- map_dbl(object@models, getAlpha)
    if (1L == length(alpha0)) {
      alpha <- alpha0 * fraction
    } else {
      alpha <- cumsum(alpha0[1:2] * c(fraction, 1 - fraction))
    }

    alpha
  })

#' @importFrom checkmate qassert
setReplaceMethod("setAlpha", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invAlpha(object, value)

    invisible(object)
  })



setMethod("getPartition", "H2ExCalibrationParam",
  function(object) {
    object@partition
  })
#' @include checkmate.R
setReplaceMethod("setPartition", "H2ExCalibrationParam",
  function(object, value) {
    assert_partition(value)
    object@partition <- value
    setDimension(object) <- length(unlist(value))

    invisible(object)
  })

setMethod("getFraction", "H2ExMarkovParam",
  function(object) {
    object@fraction
  })
#' @importFrom checkmate qassert
setReplaceMethod("setFraction", "H2ExMarkovParam",
  function(object, value) {
    qassert("N1[0,1]", value)
    object@fraction <- value

    invisible(object)
  })

setMethod("getModels", "H2ExMarkovParam",
  function(object) {
    object@models
  })
#' @importFrom purrr map_lgl map_int pmap
#' @importFrom checkmate test_class
setReplaceMethod("setModels", "H2ExMarkovParam",
  function(object, value) {
    assert_true(all(map_lgl(value, test_class, classes = getModelName(object))))
    dims <- getDimension(value)
    assert_true(dims[[1]] == sum(dims[-1]))
    partition <- pmap(
      dims[-1], cumsum(c(0, dims[2:(length(dims)-1)])), ~{
        .y + 1:.x
      })
    setPartition(object) <- partition
    object@models <- value

    invisible(object)
  })
