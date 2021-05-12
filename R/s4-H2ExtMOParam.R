#' @include s4-H2ExMOParam.R
NULL

#' H2-Extendible Marshall--Olkin calibration parameter
#'
#' Calibration parameter class for the general 2-level
#' hierarchically-etendiblew model from the Marshall--Olkin class.
#'
#' @export H2ExtMOParam
H2ExtMOParam <- setClass("H2ExtMOParam", # nolint
  contains = "H2ExMOParam")
