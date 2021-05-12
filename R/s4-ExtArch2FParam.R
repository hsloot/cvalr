#' @include s4-CalibrationParam.R s4-ExtMO2FParam.R checkmate.R
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


#' @importFrom copula getTheta
#' @importFrom checkmate qassert assert_true
setValidity("ExtArch2FParam",
  function(object) {
    qassert(object@lambda, "N1(0,)")
    qassert(object@nu, "N1")
    qassert(object@survival, "B1")
    assert_true(object@dim == object@copula@dimension)
    assert_true(object@nu == getTheta(object@copula))

    invisible(TRUE)
  })


#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,ExtArch2FParam-method
#' @aliases initialize,ExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param dim Dimension.
#' @param lambda Marginal intensity.
#' @param nu Dependence parameter (see [copula::archmCopula-class]).
#' @param rho Spearman's Rho.
#' @param tau Kendall's Tau.
#' @param survival Flag if survival copula is specified (default, except for Clayton)
#' @param family Name of the Archimedean copula family
#'   (see [copula::archmCopula-class]).
#'
#' @examples
#' ExtArch2FParam(dim = 2L, lambda = 8e-2, rho = 4e-1, family = "Gumbel")
#' ExtArch2FParam(dim = 2L, lambda = 8e-2, tau = 4e-1, family = "Amh")
#' @importFrom copula archmCopula
setMethod("initialize", "ExtArch2FParam",
  function(.Object, # nolint
      dim = 2, lambda = 0.1, nu = 1, rho = NULL, tau = NULL,
      survival = TRUE,
      family = c("Clayton", "Frank", "Amh", "Gumbel", "Joe")) {
    family <- match.arg(family)
    .Object@survival <- survival
    .Object@copula <- archmCopula(family = family)
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


#' @rdname ExtArch2FParam-class
#'
#' @export ClaytonExtArch2FParam
ClaytonExtArch2FParam <- setClass("ClaytonExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "claytonCopula"))

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,ClaytonExtArch2FParam-method
#' @aliases initialize,ClaytonExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' ClaytonExtArch2FParam(dim = 2L, lambda = 8e-2, rho = 4e-1)
#' ClaytonExtArch2FParam(dim = 2L, lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "ClaytonExtArch2FParam",
  function(.Object, ..., survival = FALSE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Clayton"))
  })


#' @rdname ExtArch2FParam-class
#'
#' @export FrankExtArch2FParam
FrankExtArch2FParam <- setClass("FrankExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "frankCopula"))

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,FrankExtArch2FParam-method
#' @aliases initialize,FrankExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' FrankExtArch2FParam(dim = 2L, lambda = 8e-2, rho = 4e-1)
#' FrankExtArch2FParam(dim = 2L, lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "FrankExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Frank"))
  })


#' @rdname ExtArch2FParam-class
#'
#' @export GumbelExtArch2FParam
GumbelExtArch2FParam <- setClass("GumbelExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "gumbelCopula"))

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,GumbelExtArch2FParam-method
#' @aliases initialize,GumbelExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' GumbelExtArch2FParam(dim = 2L, lambda = 8e-2, rho = 4e-1)
#' GumbelExtArch2FParam(dim = 2L, lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "GumbelExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Gumbel"))
  })


#' @rdname ExtArch2FParam-class
#'
#' @export AmhExtArch2FParam
AmhExtArch2FParam <- setClass("AmhExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "amhCopula"))

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,AmhExtArch2FParam-method
#' @aliases initialize,AmhExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' AmhExtArch2FParam(dim = 2L, lambda = 8e-2, rho = 4e-1)
#' AmhExtArch2FParam(dim = 2L, lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "AmhExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Amh"))
  })


#' @rdname ExtArch2FParam-class
#'
#' @export JoeExtArch2FParam
JoeExtArch2FParam <- setClass("JoeExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "joeCopula"))

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,JoeExtArch2FParam-method
#' @aliases initialize,JoeExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' JoeExtArch2FParam(dim = 2L, lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "JoeExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Joe"))
  })
