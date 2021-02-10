#' @include allClass.R allGeneric.R
NULL

setMethod("getDimension", "CalibrationParam",
  function(object) {
    object@dim
  })

setReplaceMethod("setDimension", "CalibrationParam",
  function(object, value) {
    stopifnot(1L == length(value), value %% 1 == 0, value > 0)
    object@dim <- as.integer(value)

    invisible(object)
  })


setMethod("getQMatrix", "ExMarkovParam",
  function(object) {
    object@qmatrix
  })

setReplaceMethod("setQMatrix", "ExMarkovParam",
  function(object, value) {
    stopifnot(nrow(value) == ncol(value),
      all(value[lower.tri(value)] == 0), all(value[upper.tri(value)] >= 0),
      isTRUE(all.equal(rep(0, nrow(value)), apply(value, 1, sum),
        tol = .Machine$double.eps^0.5)))

    dim <- nrow(value)-1
    setDimension(object) <- dim
    object@qmatrix <- value

    invisible(object)
  })


setMethod("getExIntensities", "ExMOParam",
  function(object) {
    object@ex_intensities
  })

setReplaceMethod("setExIntensities", "ExMOParam",
  function(object, value) {
    stopifnot(all(value >= 0), any(value > 0))
    setDimension(object) <- length(value)
    object@ex_intensities <- value
    setQMatrix(object) <- ex_intensities2qmatrix(value)

    invisible(object)
  })


setMethod("getBernsteinFunction", "ExtMOParam",
  function(object) {
    object@bf
  })

setReplaceMethod("setBernsteinFunction", "ExtMOParam",
  function(object, value) {
    stopifnot(is(value, "BernsteinFunction"))
    object@bf <- value
    setExIntensities(object) <- rmo:::bf2ex_intensities(object@dim, object@bf)

    invisible(object)
  })
#' @importFrom rmo ScaledBernsteinFunction valueOf
setReplaceMethod("setBernsteinFunction", "ExtMO2FParam",
  function(object, value) {
    stopifnot(is(value, "ScaledBernsteinFunction"),
      isTRUE(all.equal(1, valueOf(value@original, 1, 0L))))
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

setReplaceMethod("setLambda", "ExtMO2FParam",
  function(object, value) {
    stopifnot(1L == length(value), 0 < value)
    setBernsteinFunction(object) <- constructBernsteinFunction(
      object, value, object@nu)

    invisible(object)
  })
setReplaceMethod("setLambda", "ExtGaussian2FParam",
  function(object, value) {
    stopifnot(1L == length(value), 0 < value)
    object@lambda <- value

    invisible(object)
  })
setReplaceMethod("setLambda", "ExtArch2FParam",
  function(object, value) {
    stopifnot(1L == length(value), 0 < value)
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

setReplaceMethod("setNu", "ExtMO2FParam",
  function(object, value) {
    stopifnot(1L == length(value))
    setBernsteinFunction(object) <- constructBernsteinFunction(
      object, object@lambda, value)

    invisible(object)
  })
setReplaceMethod("setNu", "ExtGaussian2FParam",
  function(object, value) {
    stopifnot(1L == length(value), 0 <= value, value <= 1)
    object@nu <- value

    invisible(object)
  })
setReplaceMethod("setNu", "ExtArch2FParam",
  function(object, value) {
    stopifnot(1L == length(value))
    object@nu <- value

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
setMethod("getRho", "FrankExtArch2FParam",
  function(object) {
    copula::rho(frankCopula(object@nu))
  })

setReplaceMethod("setRho", "ExtMO2FParam",
  function(object, value) {
    stopifnot(0 <= value, value <= 1)
    setNu(object) <- invRho(object, value)

    invisible(object)
  })
setReplaceMethod("setRho", "ExtGaussian2FParam",
  function(object, value) {
    stopifnot(0 <= value, value <= 1)
    setNu(object) <- invRho(object, value)

    invisible(object)
  })
setReplaceMethod("setRho", "FrankExtArch2FParam",
  function(object, value) {
    stopifnot(0 <= value, value <= 1)
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
setMethod("getTau", "FrankExtArch2FParam",
  function(object) {
    copula::tau(frankCopula(object@nu))
  })

setReplaceMethod("setTau", "ExtMO2FParam",
  function(object, value) {
    stopifnot(0 <= value, value <= 1)
    setNu(object) <- invTau(object, value)

    invisible(object)
  })
setReplaceMethod("setTau", "ExtGaussian2FParam",
  function(object, value) {
    stopifnot(0 <= value, value <= 1)
    setNu(object) <- invTau(object, value)

    invisible(object)
  })
setReplaceMethod("setTau", "FrankExtArch2FParam",
  function(object, value) {
    stopifnot(0 <= value, value <= 1)
    setNu(object) <- invTau(object, value)

    invisible(object)
  })


#' @importFrom rmo valueOf
setMethod("getAlpha", "ExtMO2FParam",
  function(object) {
    2 - valueOf(object@bf, 2, 0L) / valueOf(object@bf, 1, 0L)
  })
setReplaceMethod("setAlpha", "ExtMO2FParam",
  function(object, value) {
    stopifnot(0 <= value, value <= 1)
    setNu(object) <- invAlpha(object, value)

    invisible(object)
  })
