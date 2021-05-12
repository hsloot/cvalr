#' @importFrom methods setClass setValidity setGeneric setMethod
#'   setReplaceMethod validObject is as new callNextMethod
#' @include checkmate.R
NULL

#' Virtual super-class for calibration parameters
#'
#' A virtual superclass for all calibration parameters for multivariate
#' (portfolio) credit models.
#'
#' @slot dim The dimension (number of portfolio items)
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
#'
#' @export
setClass("CalibrationParam", # nolint
  contains = "VIRTUAL",
  slots = c(dim = "integer"))


setGeneric("getDimension",
  function(object) {
    standardGeneric("getDimension")
  })
setGeneric("setDimension<-",
  function(object, value) {
    standardGeneric("setDimension<-")
  })


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
setValidity("CalibrationParam",
  function(object) {
    qassert(object@dim, "I1(0,)")

    invisible(TRUE)
  })


#' @describeIn CalibrationParam-class
#'    simulates the default times \eqn{(\tau_1, \ldots, \tau_d)} and returns a
#'    matrix `x` with `nrow(x) == n_sim` and `ncol(x) == dim(object)` if
#'    `dim(object) > 1L` and a vector `x` with `length(x) == n_sim` otherwise.
#'
#' @param object [CalibrationParam-class]-object.
#' ... Pass-through parameters.
#'
#' @export
setGeneric("simulate_dt",
  function(object, ...) {
    standardGeneric("simulate_dt")
  })

#' @describeIn CalibrationParam-class
#'   simulates the default counting process \eqn{L} and returns a matrix `x` with
#'   `nrow(x) == n_sim` and `ncol(x) == length(times)` if `length(times) > 1L`
#'   and a vector `x` with `length(x) == n_sim` otherwise.
#'
#' @param object Calibration parameter object.
#' @param times The times for which the process is supposed to be simulated.
#' @param ... Pass-through parameter.
#'
#' @export
setGeneric("simulate_adcp",
  function(object, times, ...) {
    standardGeneric("simulate_adcp")
  })

#' @describeIn CalibrationParam-class
#'    simulates the default times \eqn{(\tau_1, \ldots, \tau_d)} and returns a
#'    matrix `x` with `nrow(x) == n_sim` and `ncol(x) == dim(object)` if
#'    `dim(object) > 1L` and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_adcp,CalibrationParam-method
#'
#' @examples
#' parm <- ExMarkovParam(
#'  ex_qmatrix = matrix(
#'    c(-0.07647059, 0, 0, 0.05294118, -0.05, 0, 0.02352941, 0.05, 0),
#'    nrow = 3L, ncol = 3L))
#' simulate_adcp(parm, 1, n_sim = 5e1)
#' simulate_adcp(parm, seq(0, 5, by = 0.25), n_sim = 5e1)
#'
#' @importFrom stats rexp
#' @include utils.R
#' @export
setMethod("simulate_adcp", "CalibrationParam",
  function(object, times, ...) {
    tmp <- simulate_dt(object, ...)
    if (!is.matrix(tmp)) {
      tmp <- matrix(tmp, ncol = dim(object))
    }
    simplify2vector(dt2adcp(tmp, times))
  })
