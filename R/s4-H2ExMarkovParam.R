#' @include s4-H2ExCalibrationParam.R
NULL

#' H2-exchangeable Markovian calibration parameter
#'
#' Calibration parameter class for the general 2-level
#' hierarchically-exchangeable model with Markovian *default counting processes*
#' in the exchangeable sub-components.
#'
#' @slot models A list with the global and component models (of the type
#'   `ExMarkovParam-class`)
#' @slot fraction The proportion associated with the global model, see details.
#'
#' @details
#' We assume that \eqn{\tau} has the stochastic representation to be the
#' component-wise minimum of a global exchangeable Markovian-vector
#' \eqn{\tau^{(0)}} and a vector \eqn{(\tau^{(1)}, \ldots, \tau^{(J)})} with
#' independent exchangeable Markovian-vector sub-vectors \eqn{\tau^{(j)}}.
#'
#' @export H2ExMarkovParam
H2ExMarkovParam <- setClass("H2ExMarkovParam", # nolint
  contains = "H2ExCalibrationParam",
  slots = c(models = "list", fraction = "numeric"))

setGeneric("getFraction",
  function(object) {
    standardGeneric("getFraction")
  })
setGeneric("setFraction<-",
  function(object, value) {
    standardGeneric("setFraction<-")
  })

setGeneric("getModelName",
  function(object) {
    standardGeneric("getModelName")
  })

setGeneric("getModels",
  function(object) {
    standardGeneric("getModels")
  })
setGeneric("setModels<-",
  function(object, value) {
    standardGeneric("setModels<-")
  })
