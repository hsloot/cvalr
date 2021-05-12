#' @include s4-CalibrationParam.R
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
