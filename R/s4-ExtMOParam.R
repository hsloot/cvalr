#' @include s4-ExMOParam.R checkmate.R
NULL

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
#' @export ExtMOParam
ExtMOParam <- setClass("ExtMOParam", # nolint
  contains = "ExMOParam",
  slots = c(bf = "BernsteinFunction"))


setGeneric("getBernsteinFunction",
  function(object) {
    standardGeneric("getBernsteinFunction")
  })
setGeneric("setBernsteinFunction<-",
 function(object, value) {
   standardGeneric("setBernsteinFunction<-")
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


#' @importFrom checkmate assert_class
setValidity("ExtMOParam",
 function(object) {
   assert_class(object@bf, "BernsteinFunction")

   invisible(TRUE)
 })
