#' @include allClass.R allGeneric.R checkmate.R
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


setMethod("getQMatrix", "ExMarkovParam",
  function(object) {
    object@qmatrix
  })

setReplaceMethod("setQMatrix", "ExMarkovParam",
  function(object, value) {
    assert_qmatrix(value, min.rows=1, min.cols=1)

    dim <- nrow(value)-1
    setDimension(object) <- dim
    object@qmatrix <- value

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
    setQMatrix(object) <- ex_intensities2qmatrix(value)

    invisible(object)
  })


setMethod("getBernsteinFunction", "ExtMOParam",
  function(object) {
    object@bf
  })

#' @importFrom checkmate assert_class
setReplaceMethod("setBernsteinFunction", "ExtMOParam",
  function(object, value) {
    assert_class(value, "BernsteinFunction")
    object@bf <- value
    setExIntensities(object) <- rmo:::bf2ex_intensities(object@dim, object@bf)

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
#' @importFrom checkmate qassert
setReplaceMethod("setNu", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1")
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
setReplaceMethod("setRho", "FrankExtArch2FParam",
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
setMethod("getTau", "FrankExtArch2FParam",
  function(object) {
    copula::tau(frankCopula(object@nu))
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
setReplaceMethod("setTau", "FrankExtArch2FParam",
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
#' @importFrom checkmate qassert
setReplaceMethod("setAlpha", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invAlpha(object, value)

    invisible(object)
  })
