#' @include s4-H2ExCalibrationParam.R s4-H2ExtMO3FParam.R
NULL

#' Three-factor H2-extendible Gaussian calibration parameter classes
#'
#' Calibration parameter classes with three parameters for the 2-level
#' hierarchically-extendible Gaussian family with exponential margins.
#'
#' @slot lambda The marginal rate
#' @slot nu Model-specific dependence parameters for the global- and the
#'   component-models.
#'
#' @details
#' For all implemented families, the parameters `nu` can be replaced by
#' *Spearman's Rho* `rho`, *Kendall's Tau* `tau`, or the *(lower) tail
#' dependence coefficient* `alpha`.
#' For all implemented families, the possible range for `rho`, `tau`, and
#' `alpha` is between zero and one with the restriction that the global
#' parameter has to be smaller or equal than the corresponding component
#' parameter. Additionally, we support only that parameters of the same type are
#' provided, i.e. `rho`. The parameters have a one-to-one mapping to `nu`.
#'
#' @export H2ExtGaussian3FParam
H2ExtGaussian3FParam <- setClass("H2ExtGaussian3FParam", # nolint
  contains = "H2ExCalibrationParam",
  slots = c(lambda = "numeric", nu = "numeric"))
