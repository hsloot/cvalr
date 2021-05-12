#' @include s4-H2ExCalibrationParam.R
NULL

#' Three-factor H2-extendible Archimedean calibration parameter classes
#'
#' Calibration parameter classes with three parameters for the 2-level
#' hierarchically-extendible Archimedean families with exponential margins.
#'
#' @slot lambda The marginal rate
#' @slot nu Model-specific dependence parameters for the global- and the
#'   component-models.
#' @slot partition Partition of the components (only adjacent grouping allowed)
#'
#' @details
#' For all implemented families, the parameters `nu` can be replaced by
#' *Spearman's Rho* `rho` or *Kendall's Tau* `tau`.
#' For all implemented families, the possible range for `rho` and `tau` is
#' between zero and one with the restriction that the global parameter has to be
#' smaller or equal than the corresponding component parameter. Additionally, we
#' support only that parameters of the same type are provided, i.e. `rho`. The
#' parameters have a one-to-one mapping to `nu`.
#'
#' @export H2ExtArch3FParam
H2ExtArch3FParam <- setClass("H2ExtArch3FParam", # nolint
  contains = c("H2ExCalibrationParam"),
  slots = c(lambda = "numeric", nu = "numeric", family = "character",
    survival = "logical", copula = "outer_nacopula"))

#' @rdname H2ExtArch3FParam-class
#'
#' @export ClaytonH2ExtArch3FParam
ClaytonH2ExtArch3FParam <- setClass("ClaytonH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")


#' @rdname H2ExtArch3FParam-class
#'
#' @export FrankH2ExtArch3FParam
FrankH2ExtArch3FParam <- setClass("FrankH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")


#' @rdname H2ExtArch3FParam-class
#'
#' @export GumbelH2ExtArch3FParam
GumbelH2ExtArch3FParam <- setClass("GumbelH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")


#' @rdname H2ExtArch3FParam-class
#'
#' @export AmhH2ExtArch3FParam
AmhH2ExtArch3FParam <- setClass("AmhH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")


#' @rdname H2ExtArch3FParam-class
#'
#' @export JoeH2ExtArch3FParam
JoeH2ExtArch3FParam <- setClass("JoeH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")
