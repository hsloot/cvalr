#' @include s4-ExtMOParam.R checkmate.R
NULL

#' Two-factor extendible Marshall--Olkin calibration parameter classes
#'
#' Calibration parameter classes with two parameters for the extendible
#' Marshall--Olkin family.
#'
#' @slot lambda The marginal rate
#' @slot nu Model specific dependence parameter
#'
#' @details
#' For all implemented families, the parameter `nu` can be replaced by
#' *Spearman's Rho* `rho`, *Kendall' Tau* `tau` or the
#' *(lower) tail dependence coefficient* `alpha`.
#' For all implemented families, the possible range for `rho`, `tau`, and
#' `alpha` is from zero to one (boundaries might not be included) and have a
#' one-to-one mapping to the model-specific parameter `nu`.
#' The link between lower tail dependence coefficient \eqn{\alpha} and
#' Spearman's Rho and Kendall's Tau is
#' \itemize{
#'   \item \eqn{\alpha = 4 \rho / (3 + \rho)} and \eqn{\rho = 3 \alpha / (4 - \alpha)}
#'   \item \eqn{\alpha = 2 \tau / (1 + \tau)} and \eqn{\tau = \alpha / (2 - \alpha)}
#' }
#'
#' @export
setClass("ExtMO2FParam", # nolint
  contains = c("ExtMOParam", "VIRTUAL"),
  slots = c(lambda = "numeric", nu = "numeric"))


setGeneric("getLambda",
  function(object) {
    standardGeneric("getLambda")
  })
setGeneric("setLambda<-",
  function(object, value) {
    standardGeneric("setLambda<-")
  })

setGeneric("getNu",
  function(object) {
    standardGeneric("getNu")
  })
setGeneric("setNu<-",
  function(object, value) {
    standardGeneric("setNu<-")
  })

setGeneric("getRho",
  function(object) {
    standardGeneric("getRho")
  })
setGeneric("setRho<-",
  function(object, value) {
    standardGeneric("setRho<-")
  })

setGeneric("getTau",
  function(object) {
    standardGeneric("getTau")
  })
setGeneric("setTau<-",
  function(object, value) {
    standardGeneric("setTau<-")
  })

setGeneric("getAlpha",
  function(object) {
    standardGeneric("getAlpha")
  })
setGeneric("setAlpha<-",
  function(object, value) {
    standardGeneric("setAlpha<-")
  })

setGeneric("invRho",
  function(object, value) {
    standardGeneric("invRho")
  })
setGeneric("invTau",
  function(object, value) {
    standardGeneric("invTau")
  })
setGeneric("invAlpha",
  function(object, value) {
    standardGeneric("invAlpha")
  })

setGeneric("constructBernsteinFunction",
  function(object, ...) {
    standardGeneric("constructBernsteinFunction")
  })


setMethod("getLambda", "ExtMO2FParam",
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

setMethod("getNu", "ExtMO2FParam",
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

setMethod("getRho", "ExtMO2FParam",
  function(object) {
    alpha <- getAlpha(object)

    3 * alpha / (4 - alpha)
  })
#' @importFrom checkmate qassert
setReplaceMethod("setRho", "ExtMO2FParam",
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
#' @importFrom checkmate qassert
setReplaceMethod("setTau", "ExtMO2FParam",
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


#' @importFrom rmo valueOf ScaledBernsteinFunction
#' @importFrom checkmate assert qassert check_choice check_class
setValidity("ExtMO2FParam",
  function(object) {
    qassert(object@lambda, "N1(0,)")
    qassert(object@nu, "N1")
    assert(combine = "and",
      check_class(object@bf, "ScaledBernsteinFunction"),
      check_choice(object@bf@scale, object@lambda),
      check_equal(1, valueOf(object@bf@original, 1, 0L)),
      check_equal(
        object@nu, invAlpha(object, 2 - valueOf(object@bf@original, 2, 0L))))

    invisible(TRUE)
  })


#' @describeIn ExtMO2FParam-class Constructor
#' @aliases initialize,ExtMO2FParam-method
#' @aliases initialize,ExtMO2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param dim Dimension.
#' @param lambda Marginal intensity.
#' @param nu Dependence parameter.
#' @param rho Spearman's Rho.
#' @param tau Kendall's Tau.
#' @param alpha Bivariate lower tail dependence coefficient
#'
#' @examples
#' CuadrasAugeExtMO2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
#' AlphaStableExtMO2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
#' PoissonExtMO2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
#' ExponentialExtMO2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
setMethod("initialize", signature = "ExtMO2FParam", # nolint
  definition = function(.Object, # nolint
      dim, lambda, nu, rho = NULL, tau = NULL, alpha = NULL) {
    if (!missing(dim) && !missing(lambda) &&
          !(missing(nu) && missing(rho) && missing(tau) && missing(alpha))) {
      if (missing(nu)) {
        if (!is.null(rho)) {
          nu <- invRho(.Object, rho)
        } else if (!is.null(tau)) {
          nu <- invTau(.Object, tau)
        } else if (!is.null(alpha)) {
          nu <- invAlpha(.Object, alpha)
        }
      }

      setDimension(.Object) <- dim
      setBernsteinFunction(.Object) <- constructBernsteinFunction(.Object, lambda, nu)
      validObject(.Object)
    }

    invisible(.Object)
  })


#' @importFrom checkmate qassert
setMethod("invRho", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    invAlpha(object, 4 * value / (3 + value))
  })

#' @importFrom checkmate qassert
setMethod("invTau", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    invAlpha(object, 2 * value / (1 + value))
  })


#' @rdname ExtMO2FParam-class
#'
#' @section Cuadras-Augé calibration parameter class:
#' Corresponds to a Lévy subordinator which is a convex combination of
#' a pure-killing subordinator and a pure-drift subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = \nu + (1 - \nu) x}
#'   \item \eqn{\alpha = \nu}
#' }
#'
#' @export CuadrasAugeExtMO2FParam
CuadrasAugeExtMO2FParam <- setClass("CuadrasAugeExtMO2FParam", # nolint
  contains = "ExtMO2FParam")

#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions
#'   LinearBernsteinFunction ConstantBernsteinFunction
#' @importFrom checkmate assert check_choice check_class
setValidity("CuadrasAugeExtMO2FParam",
  function(object) {
    assert(combine = "and",
      check_choice(object@bf@scale, object@lambda),
      check_class(object@bf@original, "SumOfBernsteinFunctions"),
      check_class(object@bf@original@first, "LinearBernsteinFunction"),
      check_class(object@bf@original@second, "ConstantBernsteinFunction"))

    invisible(TRUE)
  })


#' @importFrom rmo SumOfBernsteinFunctions LinearBernsteinFunction
#'   ConstantBernsteinFunction
setMethod("constructBernsteinFunction", "CuadrasAugeExtMO2FParam",
  function(object, lambda, nu, ...) {
    ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = 1 - nu),
        second = ConstantBernsteinFunction(constant = nu))
    )
  })

#' @importFrom checkmate qassert
setMethod("invAlpha", "CuadrasAugeExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    value
  })



#' @rdname ExtMO2FParam-class
#'
#' @section Alpha-stable calibration parameter class:
#' Corresponds to an \eqn{\alpha}-stable subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = x^\nu}
#'   \item \eqn{\nu = \log_2(2 - \alpha)} and \eqn{\alpha = 2 - 2^\nu}
#' }
#'
#' @export AlphaStableExtMO2FParam
AlphaStableExtMO2FParam <- setClass("AlphaStableExtMO2FParam", # nolint
  contains = "ExtMO2FParam")

#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions
#'   AlphaStableBernsteinFunction
#' @importFrom checkmate assert check_choice check_class
setValidity("AlphaStableExtMO2FParam",
  function(object) {
    assert(combine = "and",
      check_choice(object@bf@scale, object@lambda),
      check_class(object@bf@original, "AlphaStableBernsteinFunction"))

    invisible(TRUE)
  })


#' @importFrom rmo SumOfBernsteinFunctions AlphaStableBernsteinFunction
setMethod("constructBernsteinFunction", "AlphaStableExtMO2FParam",
  function(object, lambda, nu, ...) {
    ScaledBernsteinFunction(
      scale = lambda,
      original = AlphaStableBernsteinFunction(alpha = nu)
    )
  })

#' @importFrom checkmate qassert
setMethod("invAlpha", "AlphaStableExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    log2(2 - value)
  })



#' @rdname ExtMO2FParam-class
#'
#' @section Poisson calibration parameter class:
#' Corresponds to a Lévy subrodinator which is a convex combination of a
#' Poisson subordinator with jump size `nu` and a pure-drift subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = \operatorname{e}^{-\nu}x + (1 - \operatorname{e}^{-x \nu})}
#'  \item \eqn{\nu = -log(1 - sqrt(\alpha))} and \eqn{\alpha = (1 - \operatorname{e}^{-\eta})}
#' }
#'
#' @export PoissonExtMO2FParam
PoissonExtMO2FParam <- setClass("PoissonExtMO2FParam", # nolint
  contains = "ExtMO2FParam")


#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions
#'   LinearBernsteinFunction PoissonBernsteinFunction
#' @importFrom checkmate assert check_choice check_class
setValidity("PoissonExtMO2FParam",
  function(object) {
    assert(combine = "and",
      check_choice(object@bf@scale, object@lambda),
      check_class(object@bf@original, "SumOfBernsteinFunctions"),
      check_class(object@bf@original@first, "LinearBernsteinFunction"),
      check_class(object@bf@original@second, "PoissonBernsteinFunction"))

      invisible(TRUE)
  })


#' @importFrom rmo SumOfBernsteinFunctions LinearBernsteinFunction
#'   PoissonBernsteinFunction
setMethod("constructBernsteinFunction", "PoissonExtMO2FParam",
  function(object, lambda, nu, ...) {
    ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = exp(-nu)),
        second = PoissonBernsteinFunction(lambda = 1, eta = nu)
      )
    )
  })

#' @importFrom checkmate qassert
setMethod("invAlpha", "PoissonExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    -log(1 - sqrt(value))
  })



#' @rdname ExtMO2FParam-class
#'
#' @section Exponential calibration parameter class:
#' Corresponds to a Lévy subordinator which is a convex combination of an
#' Exponential-jump compound Poisson process with rate `nu` and unit-intensity
#' and a pure-drift subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = (1 - 1 / (1 + \nu))x + 1 / (x + \nu)}
#'   \item \eqn{\nu = 0.5 \cdot (-3 + \sqrt{1 + 8 / \alpha})}
#'     and \eqn{\alpha = 2 / (1 + \nu) - 1 / (2 + \nu)}
#' }
#'
#' @export ExponentialExtMO2FParam
ExponentialExtMO2FParam <- setClass("ExponentialExtMO2FParam", # nolint
  contains = "ExtMO2FParam")

#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions
#'   LinearBernsteinFunction ExponentialBernsteinFunction
#' @importFrom checkmate assert check_choice check_class
setValidity("ExponentialExtMO2FParam",
  function(object) {
    assert(combine = "and",
      check_choice(object@bf@scale, object@lambda),
      check_class(object@bf@original, "SumOfBernsteinFunctions"),
      check_class(object@bf@original@first, "LinearBernsteinFunction"),
      check_class(object@bf@original@second, "ExponentialBernsteinFunction"))

      invisible(TRUE)
  })


#' @importFrom rmo SumOfBernsteinFunctions LinearBernsteinFunction
#'   ExponentialBernsteinFunction
setMethod("constructBernsteinFunction", "ExponentialExtMO2FParam",
  function(object, lambda, nu, ...) {
    ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = 1 - 1 / (1 + nu)),
        second = ExponentialBernsteinFunction(lambda = nu)
      )
    )
  })

#' @importFrom checkmate qassert
setMethod("invAlpha", "ExponentialExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    0.5 * (-3 + sqrt(1 + 8 / value))
  })
