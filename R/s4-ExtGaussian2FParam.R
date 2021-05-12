#' @include s4-CalibrationParam.R
NULL

#' Two-factor extendible Gaussian calibration parameter classes
#'
#' Calibration parameter classes with two parameters for the extendible
#' Gaussian equi-correlation family with exponential margins.
#'
#' @slot lambda The marginal rate
#' @slot nu Model specific parameter (Pearson correlation)
#'
#' @details
#' For all implemented families, the parameter `nu` can be replaced by
#' *Spearman's Rho* `rho`, *Kendall' Tau* `tau`.
#' For all implemented families, the possible range for `rho` and `tau`
#' is from zero to one (boundaries might not be included) and have a
#' one-to-one mapping to the model-specific parameter `nu`.
#' The link between the Pearson correlation coefficient and
#' Spearman's Rho and Kendall's Tau is
#' \itemize{
#'   \item \eqn{\rho = 2 \sin(\rho_S \cdot \pi / 6)} and
#'     \eqn{\rho_S = 6 / \pi \cdot \arcsin(\rho/2)}
#'   \item \eqn{\rho = \sin(\tau \cdot \pi / 2)} and
#'     \eqn{\tau = 2 / \pi \cdot \arcsin(\rho)}
#' }
#'
#' @export ExtGaussian2FParam
ExtGaussian2FParam <- setClass("ExtGaussian2FParam", # nolint
  contains = "CalibrationParam",
  slots = c("lambda", "nu"))
