#' @include s4-CalibrationParam.R s4-ExtMO2FParam.R
NULL

#' Two-factor extendible Archimedean calibration parameter classes
#'
#' Calibration parameter classes with two parameters for the extendible
#' Archimedean families with exponential margins.
#'
#' @slot lambda The marginal rate
#' @slot nu Model specific parameter
#'
#' @details
#' For all implemented families, the parameter `nu` can be replaced by
#' *Spearman's Rho* `rho`, *Kendall's Tau* `tau`.
#' For all implemented families, the possible range for `rho` and `tau`
#' is from zero to one (boundaries might not be included) and have a
#' one-to-one mapping to the model-specific parameter `nu`.
#'
#' @importFrom copula iTau iRho tau rho frankCopula iPsi
#'
#' @export ExtArch2FParam
ExtArch2FParam <- setClass("ExtArch2FParam", # nolint
  contains = "CalibrationParam",
  slots = c(lambda = "numeric", nu = "numeric", survival = "logical",
    copula = "archmCopula"))


#' @rdname ExtArch2FParam-class
#'
#' @export ClaytonExtArch2FParam
ClaytonExtArch2FParam <- setClass("ClaytonExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "claytonCopula"))


#' @rdname ExtArch2FParam-class
#'
#' @export FrankExtArch2FParam
FrankExtArch2FParam <- setClass("FrankExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "frankCopula"))


#' @rdname ExtArch2FParam-class
#'
#' @export GumbelExtArch2FParam
GumbelExtArch2FParam <- setClass("GumbelExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "gumbelCopula"))


#' @rdname ExtArch2FParam-class
#'
#' @export AmhExtArch2FParam
AmhExtArch2FParam <- setClass("AmhExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "amhCopula"))


#' @rdname ExtArch2FParam-class
#'
#' @export JoeExtArch2FParam
JoeExtArch2FParam <- setClass("JoeExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "joeCopula"))
