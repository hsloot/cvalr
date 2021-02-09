#' @importFrom methods setClass setValidity setGeneric setMethod
#'   setReplaceMethod validObject is as new callNextMethod
#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions
#'   ConstantBernsteinFunction LinearBernsteinFunction
#'   ExponentialBernsteinFunction PoissonBernsteinFunction
#'   AlphaStableBernsteinFunction

# #### All generic methods ####

setGeneric("getDimension",
  function(object) {
    standardGeneric("getDimension")
  })
setGeneric("setDimension<-",
  function(object, value) {
    standardGeneric("setDimension<-")
  })

setGeneric("getQMatrix",
  function(object) {
    standardGeneric("getQMatrix")
  })
setGeneric("setQMatrix<-",
  function(object, value) {
    standardGeneric("setQMatrix<-")
  })

setGeneric("getExIntensities",
  function(object) {
    standardGeneric("getExIntensities")
  })
setGeneric("setExIntensities<-",
  function(object, value) {
    standardGeneric("setExIntensities<-")
  })

setGeneric("getBernsteinFunction",
  function(object) {
    standardGeneric("getBernsteinFunction")
  })
setGeneric("setBernsteinFunction<-",
 function(object, value) {
   standardGeneric("setBernsteinFunction<-")
 })

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

setGeneric("iRho",
  function(object, value) {
    standardGeneric("iRho")
  })
setGeneric("iTau",
  function(object, value) {
    standardGeneric("iTau")
  })
setGeneric("iAlpha",
  function(object, value) {
    standardGeneric("iAlpha")
  })

#' Probability vector and expected value for calibration parameter
#'
#' Calculates the probability vector and expected values of the average
#' default counting process \eqn{L}.
#'
#' @param object The calibration parameter object
#' @param t Point-in-time
#'
#' @docType methods
#' @export
setGeneric("probability_vector",
  function(object, t) {
    standardGeneric("probability_vector")
  })

#' Expected loss for calibration parameter
#'

#' @rdname probability_vector
#'
#' @param g Transformation function
#' @param ... Further arguments to `g`
#'
#' @details
#' Calculates for a function \eqn{g} and the average default counting process
#' \eqn{L} the expectation
#' \deqn{
#'   \mathbb{E}[g(L_t)] .
#' }
#' For a portfolio CDS choose \eqn{g(x) = (1 - R) x} and for a CDO tranche with
#' attachment points \eqn{l < u} and choose
#' \eqn{g(x) = min{\{ \max{\{ (1 - R) x - l, 0 \}}, u - l \}}}, where \eqn{R} is
#' the recovery rate.
#'
#' @docType methods
#' @export
setGeneric("expected_value",
  function(object, t, g, ...) {
    standardGeneric("expected_value")
  })

#' @rdname probability_vector
#'
#' @param recovery_rate The recovery rate of the portfolio CDS/CDO
#' @param ... Further arguments
#'
#' @export
setGeneric("expected_pcds_loss",
  function(object, t, recovery_rate, ...) {
    standardGeneric("expected_pcds_loss")
  })

#' @rdname probability_vector
#'
#' @param lower Lower attachment point of the CDO tranche
#' @param upper Upper attachment point of the CDO tranche
#'
#' @export
setGeneric("expected_cdo_loss",
  function(object, t, recovery_rate, lower, upper, ...) {
    standardGeneric("expected_cdo_loss")
  })


# #### CalibrationParam ####

#' Virtual class `CalibrationParam` for calibration parameters
#'
#' A virtual superclass for all calibration parameters for (portfolio) credit
#' models.
#'
#' @slot dim The dimension (no of portfolio items)
#'
#' @details
#' It is assumed that the expected value of the credit derivative
#' payment streams can be represented as a function of *expected losses* and
#' other non-model dependent quantities.
#' *Expected losses* are the expected values of the average default counting
#' process \eqn{L} under a transformation \eqn{g} on times \eqn{t_1, \ldots, t_n},
#' i.e. \eqn{E[g(L_{t_1})], \ldots, E[g(L_{t_n})]}.
#'
#' A `CalibrationParam` object is supposed to be used to abstract the calculation
#' of these *expected losses* in functions which calculate the fair value of the
#' option with the asset pricing theorem.
#'
#' @seealso [ExMarkovParam-class] [ExMOParam-class] [ExtMOParam-class]
#'   [CuadrasAugeExtMO2FParam-class] [AlphaStableExtMO2FParam-class]
#'   [PoissonExtMO2FParam-class] [ExponentialExtMO2FParam-class]
setClass("CalibrationParam", # nolint
  contains = "VIRTUAL",
  slots = c(dim = "integer"))

setValidity("CalibrationParam",
  function(object) {
    stopifnot(object@dim >= 0L)

    invisible(TRUE)
  })

setMethod("getDimension", "CalibrationParam",
  function(object) {
    object@dim
  })
setReplaceMethod("setDimension", c("CalibrationParam", "integer"),
  function(object, value) {
    stopifnot(1L == length(value), value >= 1L)
    object@dim <- value

    invisible(object)
  })
setReplaceMethod("setDimension", c("CalibrationParam", "numeric"),
  function(object, value) {
    stopifnot(1L == length(value), value %% 1 == 0)
    setDimension(object) <- as.integer(value)

    invisible(object)
  })

#' @rdname probability_vector
#' @aliases expected_value,CalibrationParam
#'
#' @examples
#' expected_value(CuadrasAugeExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3,
#'   function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2))
#' expected_value(AlphaStableExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3,
#'   function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2))
#' expected_value(PoissonExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3,
#'   function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2))
#' expected_value(ExponentialExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3,
#'   function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2))
#'
#' @export
setMethod("expected_value", "CalibrationParam",
  function(object, t, g, ...) {
    mu <- sapply(0:object@dim, function(k) g(k / object@dim, ...))

    as.vector(probability_vector(object, t) %*% mu)
  })

#' @rdname probability_vector
#' @aliases expected_pcds_loss,CalibrationParam
#'
#' @examples
#' expected_pcds_loss(CuadrasAugeExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4)
#' expected_pcds_loss(AlphaStableExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4)
#' expected_pcds_loss(PoissonExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4)
#' expected_pcds_loss(ExponentialExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4)
#'
#' @export
setMethod("expected_pcds_loss", "CalibrationParam",
  function(object, t, recovery_rate, ...) {
    mu <- sapply(
      0:object@dim,
      function(k) {
        (1 - recovery_rate) * k / object@dim
      })

    as.vector(probability_vector(object, t) %*% mu)
  })

#' @rdname probability_vector
#' @aliases expected_cdo_loss,CalibrationParam
#'
#' @examples
#' expected_cdo_loss(CuadrasAugeExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4, 0.1, 0.2)
#' expected_cdo_loss(AlphaStableExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4, 0.1, 0.2)
#' expected_cdo_loss(PoissonExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4, 0.1, 0.2)
#' expected_cdo_loss(ExponentialExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4, 0.1, 0.2)
#'
#' @export
setMethod("expected_cdo_loss", "CalibrationParam",
function(object, t, recovery_rate, lower, upper, ...) {
  mu <- sapply(
    0:object@dim,
    function(k) {
      pmin(pmax((1 - recovery_rate) * k / object@dim - lower, 0), upper - lower)
    })

  as.vector(probability_vector(object, t) %*% mu)
})

# #### ExMarkovParam ####

#' Exchangeable Markovian calibration parameter
#'
#' Calibration parameter class for the general exchangeable model with a
#' Markovian *default counting process*.
#'
#' @slot qmatrix The \eqn{(d+1) \times (d+1)} Markov-generator matrix of the
#'   default counting process
#'
#' @details
#' The probability of \eqn{j > i} portfolio items being defaulted at time
#' \eqn{t > s} conditioned on \eqn{i} portfolio items being defaulted at time
#' \eqn{s} is
#'
#' \deqn{
#'   \mathbb{P}(Z_t = j \mid Z_s = i)
#'     = \delta_{i}^\top \operatorname{e}^{(t-s) Q} \delta_{j} .
#' }
#'
#' @examples
#' ExMarkovParam(
#'  qmatrix = matrix(
#'    c(-0.07647059, 0, 0, 0.05294118, -0.05, 0, 0.02352941, 0.05, 0),
#'    nrow = 3, ncol = 3))
#'
#' @export ExMarkovParam
ExMarkovParam <- setClass("ExMarkovParam", # nolint
  contains = "CalibrationParam",
  slots = c(qmatrix = "matrix"))

setValidity("ExMarkovParam",
  function(object) {
    stopifnot(object@dim+1 == nrow(object@qmatrix),
      nrow(object@qmatrix) == ncol(object@qmatrix),
      all(object@qmatrix[lower.tri(object@qmatrix)] == 0),
      all(object@qmatrix[upper.tri(object@qmatrix)] >= 0),
      isTRUE(all.equal(rep(0, nrow(object@qmatrix)), apply(object@qmatrix, 1, sum),
                tol = .Machine$double.eps^0.5)))

    invisible(TRUE)
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

setMethod("initialize", "ExMarkovParam",
  function(.Object, # nolint
      qmatrix = 0.1 * 0.5 * matrix(c(-3, 0, 0, 2, -2, 0, 1, 2, 0), nrow=3, ncol=3)) {
    setQMatrix(.Object) <- qmatrix
    validObject(.Object)

    invisible(.Object)
  })

#' @rdname probability_vector
#' @aliases probability_vector,ExMarkovParam
#'
#' @examples
#' probability_vector(CuadrasAugeExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#' probability_vector(AlphaStableExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#' probability_vector(PoissonExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#' probability_vector(ExponentialExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#'
#' @importFrom expm expm
#' @export
setMethod("probability_vector", "ExMarkovParam",
  function(object, t) {
    as.vector(expm(t * object@qmatrix)[1, ])
  })


# #### ExMOParam ####

#' Exchangeable Marshall--Olkin calibration parameter
#'
#' Calibration parameter classfor the general exchangeable model from the
#' Marshall--Olkin class.
#'
#' @slot ex_intensities The exchangeable intensities (see details)
#'
#' @details
#' The joint survival function of all portfolio items is assumed to be
#' \deqn{
#'   P(\tau > t)
#'     = \exp{(- a_0 t_{[1]} - \cdots - a_{d-1} t_{[d]})} ,
#' }
#' for \eqn{t_{[1]} \geq \cdots \geq t_{[d]}} begin the descendingly ordered
#' version of \eqn{t} and
#' \deqn{
#'   a_{i}
#'     = \sum_{l=0}^{d-i-1} \binom{d-i-1}{l} \lambda_{l+1} .
#' }
#'
#' @examples
#' ExMOParam(ex_intensities = c(0.02647059, 0.02352941))
#'
#' @export ExMOParam
ExMOParam <- setClass("ExMOParam", # nolint
  contains = "ExMarkovParam",
  slots = c(ex_intensities = "numeric"))

setValidity("ExMOParam",
  function(object) {
    stopifnot(object@dim == length(object@ex_intensities),
      all(object@ex_intensities >= 0),
      any(object@ex_intensities > 0))

    invisible(TRUE)
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

setMethod("initialize", "ExMOParam",
  definition = function(.Object, # nolint
                        ex_intensities = 0.1 * 0.5 * c(1, 1)) {
    setExIntensities(.Object) <- ex_intensities
    validObject(.Object)

    invisible(.Object)
  })


# #### ExtMOParam ####

#' Extendible Marshall--Olkin calibration parameter
#'
#' Calibration parameter class for the general extendible model from the
#' Marshall-Olkin class.
#'
#' @slot bf Bernstein function (see details)
#'
#' @details
#' The joint survival function of all portfolio items is assumed to be
#' \deqn{
#'   P(\tau > t)
#'     = \exp{(- a_{0} t_{[1]} - \cdots - a_{d-1} t_{[d]})} ,
#' }
#' for \eqn{t_{[1]} \geq \cdots \geq t_{[d]}} begin the descendingly ordered
#' version of \eqn{t} and
#' \deqn{
#'   a_{i}
#'     = \psi{(i+1)} - \psi{(i)} .
#' }
#'
#' @examples
#' ExtMOParam(
#'   dim = 2,
#'   bf = rmo::ScaledBernsteinFunction(
#'     scale = 0.05,
#'     original = rmo::SumOfBernsteinFunctions(
#'       first = rmo::ConstantBernsteinFunction(constant = 0.4),
#'       second = rmo::LinearBernsteinFunction(scale = 1 - 0.4))
#'     ))
#'
#' @export ExtMOParam
ExtMOParam <- setClass("ExtMOParam", # nolint
  contains = "ExMOParam",
  slots = c(bf = "BernsteinFunction"))

setValidity("ExtMOParam",
  function(object) {
    invisible(TRUE)
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

setMethod("initialize", "ExtMOParam",
  definition = function(.Object, # nolint
      dim = 2L,
      bf = ScaledBernsteinFunction(
        scale = 0.05,
        original = SumOfBernsteinFunctions(
          first = ConstantBernsteinFunction(constant = 0.5),
          second = LinearBernsteinFunction(scale = 0.5)
        )
      )) {
    setDimension(.Object) <- dim
    setBernsteinFunction(.Object) <- bf
    validObject(.Object)

    invisible(.Object)
  })


# #### ExtMO2FParam ####

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
#' *(lower) tail dependen coefficient* `alpha`.
#' For all implemented families, the possible range for `rho`, `tau`, and
#' `alpha` is from zero to one (boundaries might not be included) and have a
#' one-to-one mapping to the model-specific parameter `nu`.
#' The link between lower tail dependence coefficient \eqn{\alpha} and
#' Spearman's Rho and Kendall's Tau is
#' \itemize{
#'   \item \eqn{\alpha = 4 \rho / (3 + \rho)} and \eqn{\rho = 3 \alpha / (4 - \alpha)}
#'   \item \eqn{\alpha = 2 \tau / (1 + \tau)} and \eqn{\tau = \alpha / (2 - \alpha)}
#' }
setClass("ExtMO2FParam", # nolint
  contains = c("ExtMOParam", "VIRTUAL"),
  slots = c(lambda = "numeric", nu = "numeric"))

setValidity("ExtMO2FParam",
  function(object) {
    stopifnot(1L == length(object@lambda), object@lambda > 0,
      1L == length(object@nu))

    invisible(TRUE)
  })

setMethod("getLambda", "ExtMO2FParam",
  function(object) {
    object@lambda
  })
setReplaceMethod("setLambda", "ExtMO2FParam",
  function(object, value) {
    stopifnot(1L == length(value), 0 < value)
    object@lambda <- value

    invisible(object)
  })

setMethod("getNu", "ExtMO2FParam",
  function(object) {
    object@nu
  })
setReplaceMethod("setNu", "ExtMO2FParam",
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
setReplaceMethod("setRho", "ExtMO2FParam",
  function(object, value) {
    stopifnot(0 <= value, value <= 1)

    setNu(object) <- iRho(object, value)

    invisible(object)
  })

setMethod("getTau", "ExtMO2FParam",
  function(object) {
    alpha <- getAlpha(object)

    alpha / (2 - alpha)
  })
setReplaceMethod("setTau", "ExtMO2FParam",
  function(object, value) {
    stopifnot(0 <= value, value <= 1)

    setNu(object) <- iTau(object, value)

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

    setNu(object) <- iAlpha(object, value)

    invisible(object)
  })

setMethod("iRho", "ExtMO2FParam",
  function(object, value) {
    iAlpha(object, 4 * value / (3 + value))
  })
setMethod("iTau", "ExtMO2FParam", {
  function(object, value) {
    iAlpha(object, 2 * value / (1 + value))
  }
})

#' @rdname probability_vector
#'
#' @param method Choice of method (if available)
#'
#' @examples
#' expected_pcds_loss(CuadrasAugeExtMO2FParam(dim = 75, lambda = 0.05, rho = 0.4),
#'   t = 0.25, recovery_rate = 0.4)
#' expected_pcds_loss(CuadrasAugeExtMO2FParam(dim = 75, lambda = 0.05, rho = 0.4),
#'   t = 0.25, recovery_rate = 0.4, method = "fallback")
#'
#' @export
setMethod("expected_pcds_loss", "ExtMO2FParam",
  function(object, t, recovery_rate, method = c("default", "fallback"), ...) {
    method <- match.arg(method)
    if (!isTRUE("default" == method)) {
      return(callNextMethod(object, t, recovery_rate, ...))
    } else {
      return((1 - recovery_rate) * pexp(t, rate = object@lambda))
    }
  })

# #### CuadrasAugeExtMO2FParam ####

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
#' @examples
#' CuadrasAugeExtMO2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
#'
#' @export CuadrasAugeExtMO2FParam
CuadrasAugeExtMO2FParam <- setClass("CuadrasAugeExtMO2FParam", # nolint
  contains = "ExtMO2FParam")

setValidity("CuadrasAugeExtMO2FParam",
  function(object) {
    invisible(TRUE)
  })

setMethod("iAlpha", "CuadrasAugeExtMO2FParam",
  function(object, value) {
    value
  })

setMethod("initialize", signature = "CuadrasAugeExtMO2FParam",
  definition = function(.Object, # nolint
      dim = 2, lambda = 0.1, nu = 0.5, rho = NULL, tau = NULL, alpha = NULL) {
    stopifnot(!missing(nu) || !is.null(rho) || !is.null(tau) || !is.null(alpha))

    if (missing(nu)) {
      if (!is.null(rho)) {
        nu <- iRho(.Object, rho)
      } else if (!is.null(tau)) {
        nu <- iTau(.Object, tau)
      } else if (!is.null(alpha)) {
        nu <- iAlpha(.Object, alpha)
      }
    }

    setDimension(.Object) <- dim
    setLambda(.Object) <- lambda
    setNu(.Object) <- nu
    setBernsteinFunction(.Object) <- ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = 1 - nu),
        second = ConstantBernsteinFunction(constant = nu))
    )
    validObject(.Object)

    invisible(.Object)
  })


# #### AlphaStableExtMO2FParam ####

#' @rdname ExtMO2FParam-class
#'
#' @section Alpha-stable calibration parameter class:
#' Corresponds to an \eqn{\alpha}-stable subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = x^\nu}
#'   \item \eqn{\nu = \log_2(2 - \alpha)} and \eqn{\alpha = 2 - 2^\nu}
#' }
#'
#' @examples
#' AlphaStableExtMO2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
#'
#' @export AlphaStableExtMO2FParam
AlphaStableExtMO2FParam <- setClass("AlphaStableExtMO2FParam", # nolint
  contains = "ExtMO2FParam")

setValidity("AlphaStableExtMO2FParam",
  function(object) {
    invisible(TRUE)
  })

setMethod("iAlpha", "AlphaStableExtMO2FParam",
  function(object, value) {
    log2(2 - value)
  })

setMethod("initialize", signature = "AlphaStableExtMO2FParam",
  definition = function(.Object, # nolint
      dim = 2, lambda = 0.1, nu = 0.5, rho = NULL, tau = NULL, alpha = NULL) {
    stopifnot(!missing(nu) || !is.null(rho) || !is.null(tau) || !is.null(alpha))

    if (missing(nu)) {
      if (!is.null(rho)) {
        nu <- iRho(.Object, rho)
      } else if (!is.null(tau)) {
        nu <- iTau(.Object, tau)
      } else if (!is.null(alpha)) {
        nu <- iAlpha(.Object, alpha)
      }
    }

    setDimension(.Object) <- dim
    setLambda(.Object) <- lambda
    setNu(.Object) <- nu
    setBernsteinFunction(.Object) <- ScaledBernsteinFunction(
      scale = lambda,
      original = AlphaStableBernsteinFunction(alpha = nu)
    )
    validObject(.Object)

    invisible(.Object)
  })


# #### PoissonExtMO2FParam ####

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
#' @examples
#' PoissonExtMO2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
#'
#' @export PoissonExtMO2FParam
PoissonExtMO2FParam <- setClass("PoissonExtMO2FParam", # nolint
  contains = "ExtMO2FParam")

setValidity("PoissonExtMO2FParam",
  function(object) {
    invisible(TRUE)
  })

setMethod("iAlpha", "PoissonExtMO2FParam",
  function(object, value) {
    -log(1 - sqrt(value))
  })

setMethod("initialize", signature = "PoissonExtMO2FParam",
  definition = function(.Object, # nolint
      dim = 2, lambda = 0.1, nu = 0.5, rho = NULL, tau = NULL, alpha = NULL) {
    stopifnot(!missing(nu) || !is.null(rho) || !is.null(tau) || !is.null(alpha))

    if (missing(nu)) {
      if (!is.null(rho)) {
        nu <- iRho(.Object, rho)
      } else if (!is.null(tau)) {
        nu <- iTau(.Object, tau)
      } else if (!is.null(alpha)) {
        nu <- iAlpha(.Object, alpha)
      }
    }

    setDimension(.Object) <- dim
    setLambda(.Object) <- lambda
    setNu(.Object) <- nu
    setBernsteinFunction(.Object) <- ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = exp(-nu)),
        second = PoissonBernsteinFunction(lambda = 1, eta = nu)
      )
    )
    validObject(.Object)

    invisible(.Object)
  })


# #### ExponentialExtMO2FParam ####

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
#' @examples
#' ExponentialExtMO2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
#'
#' @export ExponentialExtMO2FParam
ExponentialExtMO2FParam <- setClass("ExponentialExtMO2FParam", # nolint
  contains = "ExtMO2FParam")

setValidity("ExponentialExtMO2FParam",
  function(object) {
    invisible(TRUE)
  })

setMethod("iAlpha", "ExponentialExtMO2FParam",
  function(object, value) {
    0.5 * (-3 + sqrt(1 + 8 / value))
  })

setMethod("initialize", signature = "ExponentialExtMO2FParam",
  definition = function(.Object, # nolint
      dim = 2, lambda = 0.1, nu = 0.5, rho = NULL, tau = NULL, alpha = NULL) {
    stopifnot(!missing(nu) || !is.null(rho) || !is.null(tau) || !is.null(alpha))

    if (missing(nu)) {
      if (!is.null(rho)) {
        nu <- iRho(.Object, rho)
      } else if (!is.null(tau)) {
        nu <- iTau(.Object, tau)
      } else if (!is.null(alpha)) {
        nu <- iAlpha(.Object, alpha)
      }
    }

    setDimension(.Object) <- dim
    setLambda(.Object) <- lambda
    setNu(.Object) <- nu
    setBernsteinFunction(.Object) <- ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = 1 - 1 / (1 + nu)),
        second = ExponentialBernsteinFunction(lambda = nu)
      )
    )
    validObject(.Object)

    invisible(.Object)
  })

# #### Homogeneous Gaussian ####

#' Two-factor extendible Gaussian calibration parameter classes
#'
#' Calibration parameter classes with two parameters for the extendible
#' Gaussian equi-correlation family with exponential margins.
#'
#' @slot lambda The marginal rate
#' @slot nu Model specific parameter (Pearson correlation)
#'
#' @details
#' For all implemented families, the parameter `nu` can be replaced by
#' *Spearman's Rho* `rho`, *Kendall' Tau* `tau`.
#' For all implemented families, the possible range for `rho` and `tau`
#' is from zero to one (boundaries might not be included) and have a
#' one-to-one mapping to the model-specific parameter `nu`.
#' The link between the Pearson correlation coefficient and
#' Spearman's Rho and Kendall's Tau is
#' \itemize{
#'   \item \eqn{\rho = 2 \sin(\rho_S \cdot \pi / 6)} and
#'     \eqn{\rho_S = 6 / \pi \cdot \arcsin(\rho/2)}
#'   \item \eqn{\rho = \sin(\tau \cdot \pi / 2)} and
#'     \eqn{\tau = 2 / \pi \cdot \arcsin(\rho)}
#' }
#'
#' @examples
#' ExtGaussian2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
#'
#' @export ExtGaussian2FParam
ExtGaussian2FParam <- setClass("ExtGaussian2FParam", # nolint
  contains = "CalibrationParam",
  slots = c("lambda", "nu"))

setValidity("ExtGaussian2FParam",
  function(object) {
    stopifnot(1L == length(object@lambda), object@lambda > 0,
      1L == length(object@nu), 0 <= object@nu, object@nu <= 1)

    invisible(TRUE)
  })

setMethod("getLambda", "ExtGaussian2FParam",
  function(object) {
    object@lambda
  })
setReplaceMethod("setLambda", "ExtGaussian2FParam",
  function(object, value) {
    stopifnot(1L == length(value), value > 0)
    object@lambda <- value

    invisible(object)
  })

setMethod("getNu", "ExtGaussian2FParam",
  function(object) {
    object@nu
  })
setReplaceMethod("setNu", "ExtGaussian2FParam",
  function(object, value) {
    stopifnot(1L == length(value), 0 <= value, value <= 1)
    object@nu <- value

    invisible(object)
  })

setMethod("iRho", "ExtGaussian2FParam",
  function(object, value) {
    2 * sin(value * pi / 6)
  })

setMethod("iTau", "ExtGaussian2FParam",
  function(object, value) {
    sin(value * pi / 2)
  })

setMethod("getRho", "ExtGaussian2FParam",
  function(object) {
    (6 / pi) * asin(getNu(object) / 2)
  })
setReplaceMethod("setRho", "ExtGaussian2FParam",
  function(object, value) {
    setNu(object) <- iRho(object, value)

    invisible(object)
  })

setMethod("getTau", "ExtGaussian2FParam",
  function(object) {
    (2 / pi) * asin(getNu(object))
  })
setReplaceMethod("setTau", "ExtGaussian2FParam",
  function(object, value) {
    setNu(object) <- iTau(object, value)

    invisible(object)
  })

setMethod("initialize", signature = "ExtGaussian2FParam",
  definition = function(.Object, # nolint
      dim = 2, lambda = 0.1, nu = 0.5, rho = NULL, tau = NULL) {
    stopifnot(!missing(nu) || !is.null(rho) || !is.null(tau))

    if (missing(nu)) {
      if (!is.null(rho)) {
        nu <- iRho(.Object, rho)
      } else if (!is.null(tau)) {
        nu <- iTau(.Object, tau)
      }
    }

    setDimension(.Object) <- dim
    setLambda(.Object) <- lambda
    setNu(.Object) <- nu
    validObject(.Object)

    invisible(.Object)
  })

#' @rdname probability_vector
#' @aliases probability_vector,ExtGaussian2FParam
#'
#' @examples
#' probability_vector(ExtGaussian2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#'
#' @importFrom stats integrate pexp pnorm dnorm qnorm
#' @export
setMethod("probability_vector", "ExtGaussian2FParam",
  function(object, t) {
    sapply(0:object@dim,
      function(k) {
        int_res <- integrate(
          function(x) {
            ldp <- pnorm(
              (qnorm(pexp(t, rate = object@lambda)) - sqrt(object@nu) * x) /
              (sqrt(1 - object@nu)),
              log.p=TRUE, lower.tail = TRUE
            )
            lsp <- pnorm(
              (qnorm(pexp(t, rate = object@lambda)) - sqrt(object@nu) * x) /
              (sqrt(1 - object@nu)),
              log.p=TRUE, lower.tail = FALSE
            )
            sapply(
              exp(k * ldp + (object@dim-k) * lsp) * dnorm(x),
              function(v) {
                multiply_binomial_coefficient(v, object@dim, k)
              })
          }, lower = -Inf, upper = Inf, rel.tol = .Machine$double.eps^0.5
        )

        int_res$value
      })
  })

#' @rdname probability_vector
#' @aliases expected_pcds_loss,ExtGaussian2FParam
#'
#' @examples
#' expected_pcds_loss(ExtGaussian2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   t = 0.25, recovery_rate = 0.4)
#' expected_pcds_loss(ExtGaussian2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   t = 0.25, recovery_rate = 0.4, method = "fallback")
#'
#' @export
setMethod("expected_pcds_loss", "ExtGaussian2FParam",
  function(object, t, recovery_rate, method = c("default", "fallback"), ...) {
    method <- match.arg(method)
    if (!isTRUE("default" == method)) {
      return(callNextMethod(object, t, recovery_rate, ...))
    }

    (1 - recovery_rate) * pexp(t, rate = object@lambda)
  })

#' @rdname probability_vector
#' @aliases expected_cdo_loss,ExtGaussian2FParam
#'
#' @examples
#' expected_cdo_loss(ExtGaussian2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   t = 0.25, recovery_rate = 0.4, lower = 0.1, upper = 0.2)
#' expected_cdo_loss(ExtGaussian2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   t = 0.25, recovery_rate = 0.4, lower = 0.1, upper = 0.2, method = "fallback")
#'
#' @importFrom mvtnorm pmvnorm
#' @export
setMethod("expected_cdo_loss", "ExtGaussian2FParam",
  function(
      object, t, recovery_rate, lower, upper,
      method = c("default", "fallback"), ...) {
    method <- match.arg(method)
    if (!isTRUE("default" == method)) {
      return(callNextMethod(object, t, recovery_rate, lower, upper, ...))
    }

    stopifnot(
      1L == length(recovery_rate), 0 <= recovery_rate, recovery_rate < 1,
      1L == length(lower), 1L == length(upper),
      0 <= lower, lower < upper, upper <= 1)
    corr <- matrix(c(1, rep(-sqrt(1 - object@nu), 2), 1), nrow=2, ncol=2)
    left <- pexp(t, rate = object@lambda)
    if (lower > 0) {
      left <- pmvnorm(
        lower = rep(-Inf, 2),
        upper = c(
          -qnorm(pmin(lower / (1 - recovery_rate), 1)),
          qnorm(pexp(t, rate = object@lambda))
        ),
        corr = corr
      )
    }
    right <- pmvnorm(
      lower = rep(-Inf, 2),
      upper = c(
        -qnorm(pmin(upper / (1 - recovery_rate), 1)),
        qnorm(pexp(t, rate = object@lambda))
      ),
      corr = corr
    )

    (1 - recovery_rate) * as.numeric(left - right)
  })
