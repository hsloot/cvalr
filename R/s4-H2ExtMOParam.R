#' @include s4-H2ExMOParam.R s4-ExtMO2FParam.R
NULL

#' H2-Extendible Marshall--Olkin calibration parameter
#'
#' [CalibrationParam-class] for the H2-extendible Marshall-Olkin *(average) default counting
#' process* model. Extends [H2ExMOParam-class] and related to [ExMOParam-class].
#'
#' @export H2ExtMOParam
H2ExtMOParam <- setClass("H2ExtMOParam", # nolint
  contains = "H2ExMOParam")

setMethod("getModelName", "H2ExtMOParam",
  function(object) {
    "ExtMOParam"
  })
