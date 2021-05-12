#' @include s4-H2ExCalibrationParam.R s4-H2ExtMO3FParam.R checkmate.R
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


setMethod("getLambda", "H2ExtArch3FParam",
  function(object) {
    object@lambda
  })
#' @importFrom checkmate qassert
setReplaceMethod("setLambda", "H2ExtArch3FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    object@lambda <- value

    invisible(object)
  })

setMethod("getNu", "H2ExtArch3FParam",
  function(object) {
    object@nu
  })
#' @importFrom copula onacopulaL
#' @importFrom purrr map
#' @importFrom checkmate qassert
setReplaceMethod("setNu", "H2ExtArch3FParam",
  function(object, value) {
    qassert(value, "N2")
    object@nu <- value
    object@copula <- onacopulaL(
      family = object@family,
      nacList = list(value[[1]], c(), map(object@partition, ~list(value[[2]], .))))

    invisible(object)
  })

#' @importFrom copula rho
#' @importFrom purrr map_dbl
setMethod("getRho", "H2ExtArch3FParam",
  function(object) {
    rho <- function(x) {
      copula::rho(archmCopula(family = object@family, param = x@theta, dim = 2))
    }
    c(rho(object@copula@copula), rho(object@copula@childCops[[1]]@copula))
  })

#' @importFrom copula rho
#' @importFrom purrr map_dbl
setMethod("getTau", "H2ExtArch3FParam",
  function(object) {
    tau <- function(x) {
      copula::tau(archmCopula(family = object@family, param = x@theta, dim = 2))
    }
    c(tau(object@copula@copula), tau(object@copula@childCops[[1]]@copula))
  })


#' @importFrom checkmate qassert
setValidity("H2ExtArch3FParam",
  function(object) {
    qassert(object@lambda, "N1(0,)")
    qassert(object@nu, "N2(0,)")
    qassert(object@survival, "B1")
    assert_true(object@dim == dim(object@copula))
    assert_true(
      check_equal(
        object@nu,
        c(object@copula@copula@theta, object@copula@childCops[[1]]@copula@theta)))

    invisible(TRUE)
  })


#' @describeIn H2ExtArch3FParam-class Constructor
#' @aliases initialize,H2ExtArch3FParam-method
#' @aliases initialize,H2ExtArch3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param partition Partition of the components (only adjacent grouping allowed).
#' @param lambda Marginal intensity.
#' @param nu Dependency parameter (see [copula::archmCopula-class] and
#'   [copula::nacopula-class]).
#' @param rho Spearman's Rho.
#' @param tau Kendall's Tau.
#' @param survival Flag if survival copula is specified (default, except for Clayton)
#' @param family Name of the Archimedean copula family
#'   (see [copula::archmCopula-class]).
#'
#' @examples
#' H2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, tau = c(0.3, 0.5),
#'   survival = FALSE, family = "Clayton")
#' H2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, tau = c(0.3, 0.5),
#'   survival = TRUE, family = "Gumbel")
#'
#' @importFrom copula archmCopula onacopulaL iRho iTau
#' @importFrom purrr imap map
setMethod("initialize", "H2ExtArch3FParam",
  function(.Object, # nolint
      partition = list(1L:2L, 3L:5L),
      lambda = 1e-1, nu = c(0.2, 0.3), rho = NULL, tau = NULL,
      survival = TRUE,
      family = c("Clayton", "Frank", "AMH", "Gumbel", "Joe")) {
    family <- match.arg(family)
    .Object@family <- family
    .Object@survival <- survival
    acopula <- archmCopula(family = family)
    if (!missing(partition) && !missing(lambda) &&
          (!missing(nu) || !missing(rho) || !missing(tau))) {
      if (missing(nu)) {
        if (!is.null(rho)) {
          nu <- map_dbl(rho, ~copula::iRho(acopula, .))
        } else if (!is.null(tau)) {
          nu <- map_dbl(tau, ~copula::iTau(acopula, .))
        }
      }

      dim <- length(unlist(partition))

      .Object@dim <- as.integer(dim)
      .Object@partition <- partition
      .Object@lambda <- lambda
      .Object@nu <- nu

      .Object@copula <- onacopulaL(
        family = family, list(nu[[1]], c(), map(partition, ~list(nu[[2]], .))))

      validObject(.Object)
    }

    invisible(.Object)
  })


setMethod("getModelName", "H2ExtArch3FParam",
  function(object) {
    "ExtArch2FParam"
  })

#' @importFrom copula iRho
#' @importFrom purrr map_dbl
#' @importFrom checkmate qassert
setMethod("invRho", "H2ExtArch3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    map_dbl(value, ~copula::iRho(object@copula@copula, .))
  })

#' @importFrom copula iTau
#' @importFrom checkmate qassert
setMethod("invTau", "H2ExtArch3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    copula::iTau(object@copula@copula, value)
  })


#' @describeIn H2ExtArch3FParam-class
#'    simulates the default times \eqn{(\tau_1, \ldots, \tau_d)} and returns a
#'    matrix `x` with `nrow(x) == n_sim` and `ncol(x) == dim(object)` if
#'    `dim(object) > 1L` and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_dt,H2ExtArch3FParam-method
#'
#' @inheritParams simulate_dt
#' @param n_sim Number of samples.
#'
#' @examples
#' parm <- FrankH2ExtArch3FParam(
#'   partition = list(1:2, 3:6, 7:8),
#'   lambda = 8e-2, rho = c(0.2, 0.7))
#' simulate_dt(parm, n_sim = 5e1)
#'
#' @importFrom stats qexp
#' @importFrom copula rCopula
#' @include utils.R
setMethod("simulate_dt", "H2ExtArch3FParam",
  function(object, ..., n_sim = 1e4) {
    out <- qexp(
      rCopula(n_sim, object@copula),
      rate = object@lambda, lower.tail = !object@survival
    )

    simplify2vector(out)
  })



#' @rdname H2ExtArch3FParam-class
#'
#' @export ClaytonH2ExtArch3FParam
ClaytonH2ExtArch3FParam <- setClass("ClaytonH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")


#' @describeIn H2ExtArch3FParam-class Constructor
#' @aliases initialize,ClaytonH2ExtArch3FParam-method
#' @aliases initialize,ClaytonH2ExtArch3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' ClaytonH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, tau = c(0.3, 0.5))
#' ClaytonH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, rho = c(0.3, 0.5))
setMethod("initialize", "ClaytonH2ExtArch3FParam",
  function(.Object, ..., survival = FALSE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Clayton"))
  })


setMethod("getModelName", "ClaytonH2ExtArch3FParam",
  function(object) {
    "ClaytonExtArch2FParam"
  })



#' @rdname H2ExtArch3FParam-class
#'
#' @export FrankH2ExtArch3FParam
FrankH2ExtArch3FParam <- setClass("FrankH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")


#' @describeIn H2ExtArch3FParam-class Constructor
#' @aliases initialize,FrankH2ExtArch3FParam-method
#' @aliases initialize,FrankH2ExtArch3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' FrankH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, tau = c(0.3, 0.5))
#' FrankH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, rho = c(0.3, 0.5))
setMethod("initialize", "FrankH2ExtArch3FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Frank"))
  })


setMethod("getModelName", "FrankH2ExtArch3FParam",
  function(object) {
    "FrankExtArch2FParam"
  })



#' @rdname H2ExtArch3FParam-class
#'
#' @export GumbelH2ExtArch3FParam
GumbelH2ExtArch3FParam <- setClass("GumbelH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")


#' @describeIn H2ExtArch3FParam-class Constructor
#' @aliases initialize,GumbelH2ExtArch3FParam-method
#' @aliases initialize,GumbelH2ExtArch3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' GumbelH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, tau = c(0.3, 0.5))
#' GumbelH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, rho = c(0.3, 0.5))
setMethod("initialize", "GumbelH2ExtArch3FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Gumbel"))
  })


setMethod("getModelName", "GumbelH2ExtArch3FParam",
  function(object) {
    "GumbelExtArch2FParam"
  })



#' @rdname H2ExtArch3FParam-class
#'
#' @export AmhH2ExtArch3FParam
AmhH2ExtArch3FParam <- setClass("AmhH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")


#' @describeIn H2ExtArch3FParam-class Constructor
#' @aliases initialize,AmhH2ExtArch3FParam-method
#' @aliases initialize,AmhH2ExtArch3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' AmhH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, tau = c(0.05, 0.2))
#' AmhH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, rho = c(0.05, 0.2))
setMethod("initialize", "AmhH2ExtArch3FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "AMH"))
  })


setMethod("getModelName", "AmhH2ExtArch3FParam",
  function(object) {
    "AmhExtArch2FParam"
  })



#' @rdname H2ExtArch3FParam-class
#'
#' @export JoeH2ExtArch3FParam
JoeH2ExtArch3FParam <- setClass("JoeH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")


#' @describeIn H2ExtArch3FParam-class Constructor
#' @aliases initialize,JoeH2ExtArch3FParam-method
#' @aliases initialize,JoeH2ExtArch3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' JoeH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, tau = c(0.3, 0.5))
setMethod("initialize", "JoeH2ExtArch3FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Joe"))
  })


setMethod("getModelName", "JoeH2ExtArch3FParam",
  function(object) {
    "JoeExtArch2FParam"
  })
