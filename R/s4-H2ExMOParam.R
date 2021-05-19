#' @include s4-H2ExMarkovParam.R
NULL

#' H2-Exchangeable Marshall--Olkin calibration parameter
#'
#' [CalibrationParam-class] for the H2-exchangeable Marshall-Olkin *(average) default counting
#' process* model. Extends [H2ExMarkovParam-class] and related to [ExMOParam-class].
#'
#' @export H2ExMOParam
H2ExMOParam <- setClass("H2ExMOParam", # nolint
  contains = "H2ExMarkovParam")

setMethod("getModelName", "H2ExMOParam",
  function(object) {
    "ExMOParam"
  })
