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


#' @importFrom checkmate qassert
setReplaceMethod("setDimension", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "X1(0,)")
    object@dim <- as.integer(value)
    object@copula@dimension <- as.integer(value)

    invisible(object)
  })

setMethod("getLambda", "ExtArch2FParam",
  function(object) {
    object@lambda
  })
#' @importFrom checkmate qassert
setReplaceMethod("setLambda", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    object@lambda <- value

    invisible(object)
  })

setMethod("getNu", "ExtArch2FParam",
  function(object) {
    object@nu
  })
#' @importFrom copula setTheta
#' @importFrom checkmate qassert
setReplaceMethod("setNu", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1")
    object@nu <- value
    object@copula <- setTheta(object@copula, value)

    invisible(object)
  })

#' @importFrom copula rho frankCopula
setMethod("getRho", "ExtArch2FParam",
  function(object) {
    copula::rho(object@copula)
  })
#' @importFrom checkmate qassert
setReplaceMethod("setRho", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invRho(object, value)

    invisible(object)
  })

#' @importFrom copula tau frankCopula
setMethod("getTau", "ExtArch2FParam",
  function(object) {
    copula::tau(object@copula)
  })
#' @importFrom checkmate qassert
setReplaceMethod("setTau", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invTau(object, value)

    invisible(object)
  })


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
