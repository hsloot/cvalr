#' @include s4-H2ExCalibrationParam.R s4-H2ExtMO3FParam.R checkmate.R
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


setMethod("getLambda", "H2ExtGaussian3FParam",
  function(object) {
    object@lambda
  })
#' @importFrom checkmate qassert
setReplaceMethod("setLambda", "H2ExtGaussian3FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    object@lambda <- value

    invisible(object)
  })

setMethod("getNu", "H2ExtGaussian3FParam",
  function(object) {
    object@nu
  })
#' @importFrom checkmate assert_numeric
setReplaceMethod("setNu", "H2ExtGaussian3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, sorted = TRUE)
    object@nu <- value

    invisible(object)
  })

setMethod("getRho", "H2ExtGaussian3FParam",
  function(object) {
    (6 / pi) * asin(getNu(object) / 2)
  })

setMethod("getTau", "H2ExtGaussian3FParam",
  function(object) {
    (2 / pi) * asin(getNu(object))
  })


#' @importFrom checkmate qassert
setValidity("H2ExtGaussian3FParam",
  function(object) {
    qassert(object@lambda, "N1(0,)")
    qassert(object@nu, "N2(0,)")

    invisible(TRUE)
  })


#' @importFrom purrr imap
setMethod("initialize", "H2ExtGaussian3FParam",
  function(.Object, # nolint
      partition = list(1L:2L, 3L:5L), lambda = 1e-1, nu = c(0.2, 0.3),
      rho = NULL, tau = NULL) {
    if (!missing(partition) && !missing(lambda) &&
          (!missing(nu) || !missing(rho) || !missing(tau))) {
      if (missing(nu)) {
        if (!is.null(rho)) {
          nu <- invRho(.Object, rho)
        } else if (!is.null(tau)) {
          nu <- invTau(.Object, tau)
        }
      }

      dim <- length(unlist(partition))

      .Object@dim <- dim
      .Object@partition <- partition
      .Object@lambda <- lambda
      .Object@nu <- nu

      validObject(.Object)
    }

    invisible(.Object)
  })


setMethod("getModelName", "H2ExtGaussian3FParam",
  function(object) {
    "ExtGaussian2FParam"
  })

#' @importFrom checkmate qassert
setMethod("invRho", "H2ExtGaussian3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    2 * sin(value * pi / 6)
  })

#' @importFrom checkmate qassert
setMethod("invTau", "H2ExtGaussian3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    sin(value * pi / 2)
  })
