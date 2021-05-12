#' @include s4-CalibrationParam.R s4-ExtMO2FParam.R checkmate.R
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


setMethod("getLambda", "ExtGaussian2FParam",
  function(object) {
    object@lambda
  })
#' @importFrom checkmate qassert
setReplaceMethod("setLambda", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    object@lambda <- value

    invisible(object)
  })

setMethod("getNu", "ExtGaussian2FParam",
  function(object) {
    object@nu
  })
#' @importFrom checkmate qassert
setReplaceMethod("setNu", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    object@nu <- value

    invisible(object)
  })

setMethod("getRho", "ExtGaussian2FParam",
  function(object) {
    (6 / pi) * asin(getNu(object) / 2)
  })
#' @importFrom checkmate qassert
setReplaceMethod("setRho", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invRho(object, value)

    invisible(object)
  })

setMethod("getTau", "ExtGaussian2FParam",
  function(object) {
    (2 / pi) * asin(getNu(object))
  })
#' @importFrom checkmate qassert
setReplaceMethod("setTau", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invTau(object, value)

    invisible(object)
  })


#' @importFrom checkmate qassert
setValidity("ExtGaussian2FParam",
  function(object) {
    qassert(object@lambda, "N1(0,)")
    qassert(object@nu, "N1[0,1]")

    invisible(TRUE)
  })


#' @describeIn ExtGaussian2FParam-class Constructor
#' @aliases initialize,ExtGaussian2FParam-method
#' @aliases initialize,ExtGaussian2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param dim Dimension.
#' @param lambda Marginal intensity.
#' @param nu Dependence parameter.
#' @param rho Spearman's Rho.
#' @param tau Kendall's Tau.
#'
#' @examples
#' ExtGaussian2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
setMethod("initialize", signature = "ExtGaussian2FParam",
  definition = function(.Object, # nolint
      dim, lambda, nu, rho = NULL, tau = NULL) {
    if (!missing(dim) && !missing(lambda) &&
          !(missing(nu) && missing(rho) && missing(tau))) {
      if (missing(nu)) {
        if (!is.null(rho)) {
          nu <- invRho(.Object, rho)
        } else if (!is.null(tau)) {
          nu <- invTau(.Object, tau)
        }
      }

      setDimension(.Object) <- dim
      setLambda(.Object) <- lambda
      setNu(.Object) <- nu
      validObject(.Object)
    }

    invisible(.Object)
  })


#' @importFrom checkmate qassert
setMethod("invRho", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    2 * sin(value * pi / 6)
  })

#' @importFrom checkmate qassert
setMethod("invTau", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    sin(value * pi / 2)
  })
