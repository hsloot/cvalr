#' @include s4-CalibrationParam.R checkmate.R
NULL

#' Virtual superclass for H2-exchangeable calibration parameters
#'
#' A virtual superclass for all 2-level hierarchically-exchangeable calibration
#' parameters for multivariate (portfolio) credit models.
#'
#' @slot composition Positive integerish vector for the component-composition.
#'
#' @export
setClass("H2ExCalibrationParam",
  contains = c("CalibrationParam", "VIRTUAL"),
  slots = c(composition = "integer"))


setGeneric("getComposition",
  function(object) {
    standardGeneric("getComposition")
  })
setMethod("getComposition", "H2ExCalibrationParam",
  function(object) {
    object@composition
  })

setGeneric("setComposition<-",
  function(object, value) {
    standardGeneric("setComposition<-")
  })
#' @importFrom checkmate qassert
setReplaceMethod("setComposition", "H2ExCalibrationParam",
  function(object, value) {
    qassert(value, "X+[1,)")
    qassert(sum(value), "X1[2,)")
    object@composition <- value
    setDimension(object) <- as.integer(sum(value))

    invisible(object)
  })

setGeneric("getPartition",
  function(object) {
    standardGeneric("getPartition")
  })
#' @importFrom utils head
#' @importFrom purrr map2
setMethod("getPartition", "H2ExCalibrationParam",
  function(object) {
    composition <- getComposition(object)
    map2(c(0, cumsum(head(composition, -1L))), composition, ~{
        .x + 1:.y
      })
  })

setGeneric("setPartition<-",
  function(object, value) {
    standardGeneric("setPartition<-")
  })
#' @include checkmate.R
#' @importFrom purrr map
setReplaceMethod("setPartition", "H2ExCalibrationParam",
  function(object, value) {
    assert_partition(value)
    setComposition(object) <- map_int(value, length)

    invisible(object)
  })


#' @importFrom checkmate qassert assert_choice
setValidity("H2ExCalibrationParam",
  function(object) {
    qassert(object@composition, "I+[1,)")
    assert_choice(sum(object@composition), object@dim)

    invisible(TRUE)
  })
