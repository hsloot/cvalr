#' @include s4-ExMarkovParam.R checkmate.R
NULL

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
#' @export ExMOParam
ExMOParam <- setClass("ExMOParam", # nolint
  contains = "ExMarkovParam",
  slots = c(ex_intensities = "numeric"))


setGeneric("getExIntensities",
  function(object) {
    standardGeneric("getExIntensities")
  })
setGeneric("setExIntensities<-",
  function(object, value) {
    standardGeneric("setExIntensities<-")
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


#' @importFrom checkmate qassert assert_numeric
setValidity("ExMOParam",
  function(object) {
    assert_numeric(object@ex_intensities, lower = 0, len = object@dim)
    qassert(max(object@ex_intensities), "N1(0,)")

    invisible(TRUE)
  })


#' @describeIn ExMOParam-class Constructor
#' @aliases initialize,ExMOParam-method
#' @aliases initialize,ExMOParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ex_intensities (Scaled) exchangeable intensities of the exchangeable
#'   Marshall-Olkin distribution.
#'
#' @examples
#' ExMOParam(ex_intensities = c(0.02647059, 0.02352941))
setMethod("initialize", "ExMOParam",
  definition = function(.Object, ex_intensities) { # nolint
    if (!missing(ex_intensities)) {
      setExIntensities(.Object) <- ex_intensities
      validObject(.Object)
    }

    invisible(.Object)
  })


#' @describeIn ExMOParam-class
#'    simulates the default times \eqn{(\tau_1, \ldots, \tau_d)} and returns a
#'    matrix `x` with `nrow(x) == n_sim` and `ncol(x) == dim(object)` if
#'    `dim(object) > 1L` and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_dt,ExMOParam-method
#'
#' @inheritParams simulate_dt
#' @param method Simulation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param n_sim Number of samples.
#'
#' @examples
#' parm <- ExMOParam(ex_intensities = c(0.02647059, 0.02352941))
#' simulate_dt(parm, n_sim = 5e1)
#'
#' @importFrom rmo rexmo_markovian
#' @include utils.R
#' @export
setMethod("simulate_dt", "ExMOParam",
  function(object, ...,
      method = c("default", "ExMOParam", "ExMarkovParam"), n_sim = 1e4) {
    method <- match.arg(method)
    if (isTRUE("default" == method)) {
      method <- "ExMOParam"
    }

    if (isTRUE("ExMOParam" == method)) {
      out <- rexmo_markovian(n_sim, object@dim, object@ex_intensities)
      out <- simplify2vector(out)
    } else {
      out <- callNextMethod(object, ..., n_sim = n_sim)
    }

    out
  })
