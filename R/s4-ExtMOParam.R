#' @include s4-ExMOParam.R checkmate.R
NULL

#' Extendible Marshall--Olkin calibration parameter
#'
#' @description
#' [CalibrationParam-class]-class for the extendible Marshall-Olkin model for
#' the *average default counting process*.
#'
#' @slot bf The Bernstein function of the extendible Marshall-Olkin distribution
#' (see [rmo::BernsteinFunction-class]).
#'
#' @details
#' The model is defined by the assumption that the multivariate default times
#' \eqn{\tau = (\tau_1, \ldots, \tau_d)} are extendible Marshall-Olkin.
#' The joint survival function of all portfolio items is
#' \deqn{
#'   P(\tau > t)
#'     = \exp{(- a_0 t_{[1]} - \cdots - a_{d-1} t_{[d]})} ,
#' }
#' for \eqn{t_{[1]} \geq \cdots \geq t_{[d]}} begin the descending version of
#' \eqn{t} and
#' \deqn{
#'   a_{i}
#'     = \sum_{l=0}^{d-i-1} \binom{d-i-1}{l} \lambda_{l+1} .
#' }
#' The parameter are implicitly defined by a +Bernstein function* \eqn{\psi}
#' (which is provided to the constructor):
#' \deqn{
#'   a_{i}
#'     = \psi{(i+1)} - \psi{(i)} .
#' }
#'
#' @export ExtMOParam
ExtMOParam <- setClass("ExtMOParam", # nolint
  contains = "ExMOParam",
  slots = c(bf = "BernsteinFunction"))


setGeneric("getBernsteinFunction",
  function(object) {
    standardGeneric("getBernsteinFunction")
  })
setMethod("getBernsteinFunction", "ExtMOParam",
 function(object) {
   object@bf
 })

setGeneric("setBernsteinFunction<-",
 function(object, value) {
   standardGeneric("setBernsteinFunction<-")
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


#' @importFrom checkmate assert_class
setValidity("ExtMOParam",
 function(object) {
   assert_class(object@bf, "BernsteinFunction")

   invisible(TRUE)
 })


#' @describeIn ExtMOParam-class Constructor
#' @aliases initialize,ExtMOParam-method
#' @aliases initialize,ExtMOParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param dim Dimension.
#' @param bf Bernstein function.
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
setMethod("initialize", "ExtMOParam",
 definition = function(.Object, dim, bf) { # nolint
   if (!missing(dim) && !missing(bf)) {
     setDimension(.Object) <- dim
     setBernsteinFunction(.Object) <- bf
     validObject(.Object)
   }

   invisible(.Object)
 })


#' @describeIn ExtMOParam-class Display the object.
#' @aliases show,ExtMOParam-method
#'
#' @inheritParams methods::show
#'
#' @export
setMethod("show", "ExtMOParam",
 function(object) {
   cat("An object of class \"ExtMOParam\"\n")
   cat(sprintf("Dimension: %i\n", getDimension(object)))
   cat("Bernstein function:\n")
   print(getBernsteinFunction(object))

   invisible(NULL)
 })
