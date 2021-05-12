#' @include s4-ExMarkovParam.R
NULL

#' Exchangeable Marshall--Olkin calibration parameter
#'
#' Calibration parameter classfor the general exchangeable model from the
#' Marshall--Olkin class.
#'
#' @slot ex_intensities The exchangeable intensities (see details)
#'
#' @details
#' The joint survival function of all portfolio items is assumed to be
#' \deqn{
#'   P(\tau > t)
#'     = \exp{(- a_0 t_{[1]} - \cdots - a_{d-1} t_{[d]})} ,
#' }
#' for \eqn{t_{[1]} \geq \cdots \geq t_{[d]}} begin the descendingly ordered
#' version of \eqn{t} and
#' \deqn{
#'   a_{i}
#'     = \sum_{l=0}^{d-i-1} \binom{d-i-1}{l} \lambda_{l+1} .
#' }
#'
#' @export ExMOParam
ExMOParam <- setClass("ExMOParam", # nolint
  contains = "ExMarkovParam",
  slots = c(ex_intensities = "numeric"))
