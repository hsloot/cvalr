#' @include s4-H2ExMarkovParam.R
NULL

#' H2-Exchangeable Marshall--Olkin calibration parameter
#'
#' Calibration parameter class for the general 2-level
#' hierarchically-exchangeable model from the Marshall--Olkin class.
#'
#' @export H2ExMOParam
H2ExMOParam <- setClass("H2ExMOParam", # nolint
  contains = "H2ExMarkovParam")
