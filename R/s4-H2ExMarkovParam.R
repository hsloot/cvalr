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


setMethod("getFraction", "H2ExMarkovParam",
  function(object) {
    object@fraction
  })
#' @importFrom checkmate qassert
setReplaceMethod("setFraction", "H2ExMarkovParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    object@fraction <- value

    invisible(object)
  })

setMethod("getModels", "H2ExMarkovParam",
  function(object) {
    object@models
  })
#' @importFrom purrr map_lgl map_int map2
#' @importFrom checkmate test_class
setReplaceMethod("setModels", "H2ExMarkovParam",
  function(object, value) {
    assert_true(all(map_lgl(value, test_class, classes = getModelName(object))))
    dims <- map_int(value, getDimension)
    assert_true(dims[[1]] == sum(dims[-1]))
    partition <- map2(
      dims[-1], cumsum(c(0, dims[2:(length(dims)-1)])), ~{
        .y + 1:.x
      })
    setPartition(object) <- partition
    object@models <- value

    invisible(object)
  })
