#' @importFrom methods setClass setValidity setGeneric setMethod
#'   setReplaceMethod validObject is as new callNextMethod
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
