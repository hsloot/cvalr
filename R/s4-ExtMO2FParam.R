#' @include s4-ExtMOParam.R checkmate.R
NULL

#' Two-factor extendible Marshall-Olkin calibration parameters
#'
#' @description
#' [CalibrationParam-class]-class with two parameters for the extendible
#' Marshall-Olkin model for the *(average) default counting process*.
#' Extends [ExtMOParam-class].
#'
#' @slot lambda A non-negative number for the marginal rate.
#' @slot nu A numeric number for the model specific dependence parameter (range
#'   depends on specific model, use `rho`, `tau`, or `alpha` to set dependence
#'   parameter).
#'
#' @details
#' The model is defined by the assumption that the multivariate default times
#' \eqn{\tau = (\tau_1, \ldots, \tau_d)} are extendible Marshall-Olkin, see
#' [ExtMOParam-class] for the details. This class provides an interface for
#' easy-to-use, 2-factor families for this model.
#' For all implemented families, the marginal rate can be specified by `lambda`
#' and the dependence can be specified by the internal parameter `nu`.
#' For all implemented families, the (internal) dependence parameter `nu` has a
#' one-to-one relationship, and can be replaced by, *Spearman's Rho* `rho`,
#' *Kendall' Tau* `tau` or the *(lower) tail dependence coefficient* `alpha`.
#' The possible range for `rho`, `tau`, and `alpha` is from zero to one
#' (boundaries might not be included).
#' The link between lower tail-dependence coefficient \eqn{\alpha} and
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


setGeneric("constructBernsteinFunction",
  function(object, ...) {
    standardGeneric("constructBernsteinFunction")
  })

#' @importFrom rmo ScaledBernsteinFunction valueOf
#' @importFrom checkmate assert check_class
setReplaceMethod("setBernsteinFunction", "ExtMO2FParam",
  function(object, value) {
    assert(combine = "and",
      check_class(value, "ScaledBernsteinFunction"),
      check_equal(1, valueOf(value@original, 1, 0L)))

    callNextMethod()
  })

setGeneric("getLambda",
  function(object) {
    standardGeneric("getLambda")
  })
setMethod("getLambda", "ExtMO2FParam",
  function(object) {
    object@lambda
  })

setGeneric("setLambda<-",
  function(object, value) {
    standardGeneric("setLambda<-")
  })
#' @importFrom checkmate qassert qtest
setReplaceMethod("setLambda", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    object@lambda <- value
    nu <- getNu(object)
    if (qtest(nu, "N1")) {
      setBernsteinFunction(object) <- constructBernsteinFunction(object, value, nu)
    }

    invisible(object)
  })

setGeneric("getNu",
  function(object) {
    standardGeneric("getNu")
  })
setMethod("getNu", "ExtMO2FParam",
  function(object) {
    object@nu
  })

setGeneric("setNu<-",
  function(object, value) {
    standardGeneric("setNu<-")
  })
#' @importFrom checkmate qassert qtest
setReplaceMethod("setNu", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1")
    object@nu <- value
    lambda <- getLambda(object)
    if (qtest(lambda, "N1(0,)")) {
      setBernsteinFunction(object) <- constructBernsteinFunction(object, lambda, value)
    }

    invisible(object)
  })

setGeneric("calcAlpha2Nu",
  function(object, value) {
    standardGeneric("calcAlpha2Nu")
  })

setGeneric("calcNu2Alpha",
  function(object, value) {
    standardGeneric("calcNu2Alpha")
  })

#' @importFrom checkmate qassert
setGeneric("invAlpha",
  function(object, value) {
    standardGeneric("invAlpha")
  })
#' @importFrom checkmate qassert
setMethod("invAlpha", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    calcAlpha2Nu(object, value)
  })

setGeneric("getAlpha",
  function(object) {
    standardGeneric("getAlpha")
  })
#' @importFrom rmo valueOf
setMethod("getAlpha", "ExtMO2FParam",
  function(object) {
    getNu(object) %>%
      calcNu2Alpha(object, .)
  })

setGeneric("setAlpha<-",
  function(object, value) {
    standardGeneric("setAlpha<-")
  })
#' @importFrom checkmate qassert
setReplaceMethod("setAlpha", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invAlpha(object, value)

    invisible(object)
  })

setGeneric("calcRho2Alpha",
  function(object, value) {
    standardGeneric("calcRho2Alpha")
  })
#' @importFrom checkmate qassert
setMethod("calcRho2Alpha", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    4 * value / (3 + value)
  })

setGeneric("calcRho2Nu",
  function(object, value) {
    standardGeneric("calcRho2Nu")
  })
#' @importFrom checkmate qassert
setMethod("calcRho2Nu", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    calcRho2Alpha(object, value) %>%
      calcAlpha2Nu(object, .)
  })

setGeneric("invRho",
  function(object, value) {
    standardGeneric("invRho")
  })
#' @importFrom checkmate qassert
setMethod("invRho", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    calcRho2Nu(object, value)
  })

setGeneric("calcTau2Alpha",
  function(object, value) {
    standardGeneric("calcTau2Alpha")
  })
#' @importFrom checkmate qassert
setMethod("calcTau2Alpha", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    2 * value / (1 + value)
  })

setGeneric("calcTau2Nu",
  function(object, value) {
    standardGeneric("calcTau2Nu")
  })
#' @importFrom checkmate qassert
setMethod("calcTau2Nu", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    calcTau2Alpha(object, value) %>%
      calcAlpha2Nu(object, .)
  })

setGeneric("invTau",
  function(object, value) {
    standardGeneric("invTau")
  })
#' @importFrom checkmate qassert
setMethod("invTau", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    calcTau2Nu(object, value)
  })

setGeneric("calcAlpha2Rho",
  function(object, value) {
    standardGeneric("calcAlpha2Rho")
  })
#' @importFrom checkmate qassert
setMethod("calcAlpha2Rho", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    3 * value / (4 - value)
  })

setGeneric("getRho",
  function(object) {
    standardGeneric("getRho")
  })
#' @importFrom checkmate qassert
setMethod("getRho", "ExtMO2FParam",
  function(object) {
    getAlpha(object) %>%
      calcAlpha2Rho(object, .)
  })

setGeneric("setRho<-",
  function(object, value) {
    standardGeneric("setRho<-")
  })
#' @importFrom checkmate qassert
setReplaceMethod("setRho", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invRho(object, value)

    invisible(object)
  })

setGeneric("calcAlpha2Tau",
  function(object, value) {
    standardGeneric("calcAlpha2Tau")
  })
#' @importFrom checkmate qassert
setMethod("calcAlpha2Tau", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    value / (2 - value)
  })

setGeneric("getTau",
  function(object) {
    standardGeneric("getTau")
  })
setMethod("getTau", "ExtMO2FParam",
  function(object) {
    getAlpha(object) %>%
      calcAlpha2Tau(object, .)
  })

setGeneric("setTau<-",
  function(object, value) {
    standardGeneric("setTau<-")
  })
#' @importFrom checkmate qassert
setReplaceMethod("setTau", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invTau(object, value)

    invisible(object)
  })



#' @importFrom rmo valueOf ScaledBernsteinFunction
#' @importFrom checkmate assert qassert check_choice check_class
#' @include checkmate.R
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
#' @param nu (Internal) bivariate dependence parameter.
#' @param rho Bivariate Spearman's Rho.
#' @param tau Bivariate Kendall's Tau.
#' @param alpha Bivariate lower tail-dependence coefficient.
#'
#' @examples
#' CuadrasAugeExtMO2FParam(dim = 5L, lambda = 8e-2, rho = 4e-1)
#' AlphaStableExtMO2FParam(dim = 5L, lambda = 8e-2, rho = 4e-1)
#' PoissonExtMO2FParam(dim = 5L, lambda = 8e-2, tau = 4e-1)
#' ExponentialExtMO2FParam(dim = 5L, lambda = 8e-2, alpha = 4e-1)
setMethod("initialize", signature = "ExtMO2FParam", # nolint
  definition = function(.Object, # nolint
      dim, lambda, nu, rho, tau, alpha) {
    if (!missing(dim) && !missing(lambda) &&
          !(missing(nu) && missing(rho) && missing(tau) && missing(alpha))) {
      setDimension(.Object) <- dim
      setLambda(.Object) <- lambda
      if (!missing(nu)) {
        setNu(.Object) <- nu
      } else if (!missing(rho)) {
        setRho(.Object) <- rho
      } else if (!missing(tau)) {
        setTau(.Object) <- tau
      } else if (!missing(alpha)) {
        setAlpha(.Object) <- alpha
      }

      validObject(.Object)
    }

    invisible(.Object)
  })



#' @describeIn ExtMO2FParam-class
#'   calculates the *expected value* for the *portfolio CDS loss* based on the
#'   *average default count process* for given timepoints and returns a vector
#'   `x` with `length(x) == length(times)`.
#' @aliases expected_pcds_loss,ExtMO2FParam-method
#'
#' @inheritParams expected_pcds_loss
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#'
#' @section Expected portfolio CDS loss:
#' The *expected portfolio CDS loss* for *recovery rate* \eqn{R} is calculated
#' using that
#' \deqn{
#'   \mathbb{E}[g(L_t)]
#'     = (1 - R) \cdot F(t)
#' }
#' with \eqn{g(x) = (1 - R) \cdot x} and \eqn{F} being the Exponential
#' distribution function for rate \eqn{\lambda}.
#'
#' @examples
#' expected_pcds_loss(CuadrasAugeExtMO2FParam(dim = 75L, lambda = 0.05, rho = 0.4),
#'   times = 0.25, recovery_rate = 0.4)
#' expected_pcds_loss(CuadrasAugeExtMO2FParam(dim = 75L, lambda = 0.05, rho = 0.4),
#'   times = seq(0, 5, by = 0.25), recovery_rate = 0.4)
#' expected_pcds_loss(CuadrasAugeExtMO2FParam(dim = 75L, lambda = 0.05, rho = 0.4),
#'   times = seq(0, 5, by = 0.25), recovery_rate = 0.4, method = "CalibrationParam")
#' expected_pcds_loss(CuadrasAugeExtMO2FParam(dim = 75L, lambda = 0.05, rho = 0.4),
#'   times = seq(0, 5, by = 0.25), recovery_rate = 0.4, method = "CalibrationParam",
#'   pd_args = list(method = "CalibrationParam", seed = 1623, sim_args = list(n_sim = 1e2L)))
#' expected_pcds_loss(CuadrasAugeExtMO2FParam(dim = 75L, lambda = 0.05, rho = 0.4),
#'   times = seq(0, 5, by = 0.25), recovery_rate = 0.4, method = "CalibrationParam",
#'   pd_args = list(method = "CalibrationParam", seed = 1623,
#'   sim_args = list(method = "ExMarkovParam", n_sim = 1e2L)))
#'
#' @importFrom stats pexp
#' @importFrom checkmate qassert
#' @export
setMethod("expected_pcds_loss", "ExtMO2FParam",
  function(object, times, recovery_rate, ...,
      method = c("default", "ExtMO2FParam", "CalibrationParam")) {
    method <- match.arg(method)
    if (isTRUE("default" == method || "ExtMO2FParam" == method)) {
      qassert(times, "N+[0,)")
      qassert(recovery_rate, "N1[0,1]")
      out <- (1 - recovery_rate) * pexp(times, rate = getLambda(object))
    } else {
      out <- callNextMethod(object, times, recovery_rate, ..., method = method)
    }

    out
  })


#' @describeIn ExtMOParam-class Display the object.
#' @aliases show,ExtMO2FParam-method
#'
#' @inheritParams methods::show
#'
#' @export
setMethod("show", "ExtMO2FParam",
 function(object) {
   cat(sprintf("An object of class %s\n", classLabel(class(object))))
   cat(sprintf("Dimension: %i\n", getDimension(object)))
   cat("Parameter:\n")
   cat(sprintf("* %s: %s\n", "Lambda", format(getLambda(object))))
   cat(sprintf("* %s: %s\n", "Rho", format(getRho(object))))
   cat(sprintf("* %s: %s\n", "Tau", format(getTau(object))))
   cat(sprintf("* %s: %s\n", "Alpha", format(getAlpha(object))))
   cat("Internal parameter:\n")
   cat(sprintf("* %s: %s\n", "Nu", format(getNu(object))))
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

#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions LinearBernsteinFunction
#'   ConstantBernsteinFunction
#' @importFrom checkmate assert check_choice check_class qassert
setValidity("CuadrasAugeExtMO2FParam",
  function(object) {
    assert(combine = "and",
      check_choice(object@bf@scale, object@lambda),
      check_class(object@bf@original, "SumOfBernsteinFunctions"),
      check_class(object@bf@original@first, "LinearBernsteinFunction"),
      check_class(object@bf@original@second, "ConstantBernsteinFunction"))
    qassert(object@nu, "N1[0,1]")

    invisible(TRUE)
  })


#' @importFrom rmo SumOfBernsteinFunctions LinearBernsteinFunction ConstantBernsteinFunction
#' @importFrom checkmate qassert
setMethod("constructBernsteinFunction", "CuadrasAugeExtMO2FParam",
  function(object, lambda, nu, ...) {
    qassert(lambda, "N1(0,)")
    qassert(nu, "N1[0,1]")
    ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = 1 - nu),
        second = ConstantBernsteinFunction(constant = nu))
    )
  })

#' @importFrom checkmate qassert
setMethod("calcAlpha2Nu", "CuadrasAugeExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    value
  })

#' @importFrom checkmate qassert
setMethod("calcNu2Alpha", "CuadrasAugeExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    value
  })

#' @importFrom checkmate qassert
setReplaceMethod("setNu", "CuadrasAugeExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    callNextMethod()
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

#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions AlphaStableBernsteinFunction
#' @importFrom checkmate assert check_choice check_class qassert
setValidity("AlphaStableExtMO2FParam",
  function(object) {
    assert(combine = "and",
      check_choice(object@bf@scale, object@lambda),
      check_class(object@bf@original, "AlphaStableBernsteinFunction"))
    qassert(object@nu, "N1(0,1)")

    invisible(TRUE)
  })


#' @importFrom checkmate qassert
#' @importFrom rmo SumOfBernsteinFunctions AlphaStableBernsteinFunction
setMethod("constructBernsteinFunction", "AlphaStableExtMO2FParam",
  function(object, lambda, nu) {
    qassert(lambda, "N1(0,)")
    qassert(nu, "N1(0,1)")
    ScaledBernsteinFunction(
      scale = lambda,
      original = AlphaStableBernsteinFunction(alpha = nu)
    )
  })

#' @importFrom checkmate qassert
setMethod("calcAlpha2Nu", "AlphaStableExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,1)")
    log2(2 - value)
  })

#' @importFrom checkmate qassert
setMethod("calcNu2Alpha", "AlphaStableExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,1)")
    2 - 2 ^ value
  })

setReplaceMethod("setNu", "AlphaStableExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,1)")
    callNextMethod()
  })



#' @rdname ExtMO2FParam-class
#'
#' @section Poisson calibration parameter class:
#' Corresponds to a Lévy subrodinator which is a convex combination of a
#' Poisson subordinator with jump size `nu` and a pure-drift subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = \operatorname{e}^{-\nu}x + (1 - \operatorname{e}^{-x \nu})}
#'  \item \eqn{\nu = -log(1 - sqrt(\alpha))} and \eqn{\alpha = (1 - \operatorname{e}^{-\nu})^2}
#' }
#'
#' @export PoissonExtMO2FParam
PoissonExtMO2FParam <- setClass("PoissonExtMO2FParam", # nolint
  contains = "ExtMO2FParam")


#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions LinearBernsteinFunction
#'   PoissonBernsteinFunction
#' @importFrom checkmate assert check_choice check_class qassert qassert
setValidity("PoissonExtMO2FParam",
  function(object) {
    assert(combine = "and",
      check_choice(object@bf@scale, object@lambda),
      check_class(object@bf@original, "SumOfBernsteinFunctions"),
      check_class(object@bf@original@first, "LinearBernsteinFunction"),
      check_class(object@bf@original@second, "PoissonBernsteinFunction"))
    qassert(object@nu, "N1[0,)")

      invisible(TRUE)
  })


#' @importFrom rmo SumOfBernsteinFunctions LinearBernsteinFunction PoissonBernsteinFunction
#' @importFrom checkmate qassert
setMethod("constructBernsteinFunction", "PoissonExtMO2FParam",
  function(object, lambda, nu) {
    qassert(lambda, "N1(0,)")
    qassert(nu, "N1[0,)")
    ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = exp(-nu)),
        second = PoissonBernsteinFunction(lambda = 1, eta = nu)
      )
    )
  })

#' @importFrom checkmate qassert
setMethod("calcAlpha2Nu", "PoissonExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1)")
    -log(1 - sqrt(value))
  })

#' @importFrom checkmate qassert
setMethod("calcNu2Alpha", "PoissonExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,)")
    (1 - exp(-value))^2
  })

setReplaceMethod("setNu", "PoissonExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,)")
    callNextMethod()
  })



#' @rdname ExtMO2FParam-class
#'
#' @section Exponential calibration parameter class:
#' Corresponds to a Lévy subordinator which is a convex combination of an
#' Exponential-jump compound Poisson process with rate `nu` and unit-intensity
#' and a pure-drift subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = (1 - 1 / (1 + \nu))x + x   / (x + \nu)}
#'   \item \eqn{\nu = 0.5 \cdot (-3 + \sqrt{1 + 8 / \alpha})}
#'     and \eqn{\alpha = 2 / (1 + \nu) - 2 / (2 + \nu)}
#' }
#'
#' @export ExponentialExtMO2FParam
ExponentialExtMO2FParam <- setClass("ExponentialExtMO2FParam", # nolint
  contains = "ExtMO2FParam")

#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions LinearBernsteinFunction
#'   ExponentialBernsteinFunction
#' @importFrom checkmate assert check_choice check_class qassert
setValidity("ExponentialExtMO2FParam",
  function(object) {
    assert(combine = "and",
      check_choice(object@bf@scale, object@lambda),
      check_class(object@bf@original, "SumOfBernsteinFunctions"),
      check_class(object@bf@original@first, "LinearBernsteinFunction"),
      check_class(object@bf@original@second, "ExponentialBernsteinFunction"))
    qassert(object@nu, "N1(0,)")

    invisible(TRUE)
  })


#' @importFrom rmo SumOfBernsteinFunctions LinearBernsteinFunction ExponentialBernsteinFunction
#' @importFrom checkmate qassert
setMethod("constructBernsteinFunction", "ExponentialExtMO2FParam",
  function(object, lambda, nu, ...) {
    qassert(lambda, "N1(0,)")
    qassert(nu, "N1(0,)")
    ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = 1 - 1 / (1 + nu)),
        second = ExponentialBernsteinFunction(lambda = nu)
      )
    )
  })

#' @importFrom checkmate qassert
setMethod("calcAlpha2Nu", "ExponentialExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,1)")
    0.5 * (-3 + sqrt(1 + 8 / value))
  })

#' @importFrom checkmate qassert
setMethod("calcNu2Alpha", "ExponentialExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    2 * (1 / (1 + value) - 1 / (2 + value))
  })

setReplaceMethod("setNu", "ExponentialExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    callNextMethod()
  })
