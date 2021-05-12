#' @include s4-CalibrationParam.R checkmate.R
NULL

#' Virtual superclass for H2-exchangeable calibration parameters
#'
#' A virtual superclass for all 2-level hierarchically-exchangeable calibration
#' parameters for multivariate (portfolio) credit models.
#'
#' @slot partition Partition of the components (only adjacent grouping allowed)
#'
#' @export
setClass("H2ExCalibrationParam",
  contains = c("CalibrationParam", "VIRTUAL"),
  slots = c(partition = "list"))


setGeneric("getPartition",
  function(object) {
    standardGeneric("getPartition")
  })
setGeneric("setPartition<-",
  function(object, value) {
    standardGeneric("setPartition<-")
  })


setMethod("getPartition", "H2ExCalibrationParam",
  function(object) {
    object@partition
  })
#' @include checkmate.R
#' @importFrom purrr map
setReplaceMethod("setPartition", "H2ExCalibrationParam",
  function(object, value) {
    assert_partition(value)
    object@partition <- map(value, as.integer)
    setDimension(object) <- length(unlist(value))

    invisible(object)
  })


#' @importFrom purrr map_lgl
#' @importFrom checkmate qassert qtest assert_true
setValidity("H2ExCalibrationParam",
  function(object) {
    assert_partition(object@partition)

    invisible(TRUE)
  })
