#' @include s4-H2ExtMOParam.R checkmate.R
NULL

#' Two-factor H2-extendible Marshall--Olkin calibration parameter classes
#'
#' Calibration parameter class with three parameters for the general 2-level
#' hierarchically-etendiblew model from the Marshall--Olkin class.
#'
#' @slot lambda The marginal rate
#' @slot nu Model-specific dependence parameters for the global- and the
#'   component-models.
#'
#' @details
#' For all implemented families, the parameters `nu`can be replaced by
#' *Spearman's Rho* `rho`, *Kendall's Tau* `tau`, or the *(lower) tail
#' dependence coefficient* `alpha`.
#' For all implemented families, the possible range for `rho`, `tau`, and
#' `alpha`, is between zero and one with the restriction that the global
#' parameter has to be smaller or equal than the corresponding component
#' parameter. Additionally, we support only that parameters of the same type are
#' provided, i.e. `rho`. The parameters have a one-to-one mapping to `nu`.
#'
#' @export
setClass("H2ExtMO3FParam", # nolint
  contains = c("H2ExtMOParam", "VIRTUAL"),
  slots = c(lambda = "numeric", nu = "numeric"))


setMethod("getLambda", "H2ExtMO3FParam",
  function(object) {
    object@lambda
  })
  #' @include checkmate.R
  #' @importFrom purrr map
  #' @importFrom checkmate qassert
setReplaceMethod("setLambda", "H2ExtMO3FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    assert_partition(object@partition)
    object@lambda <- value
    if (0L == length(object@models)) {
      models <- map(object@partition, ~{
        new(getModelName(object), length(.x), value, alpha = 0.5)
      })
      models <- c(list(new(getModelName(object), getDimension(object), value, alpha = 0.5)), models)
      setModels(object) <- models
    } else {
      object@models <- map(object@models, ~{
        setLambda(.x) <- value
        .x
      })
    }

    invisible(object)
  })

setMethod("getNu", "H2ExtMO3FParam",
  function(object) {
    object@nu
  })
#' @importFrom purrr map imap
#' @importFrom checkmate qassert
setReplaceMethod("setNu", "H2ExtMO3FParam",
  function(object, value) {
    qassert(value, "N2")
    assert_partition(object@partition)
    object@nu <- value
    if (0L == length(object@models)) {
      models <- map(object@partition, ~{
        new(getModelName(object), length(.x), 1, value[[2]])
      })
      models <- c(list(new(getModelName(object), getDimension(object), 1, value[[1]])), models)
      setModels(object) <- models
    } else {
      object@models <- imap(object@models, ~{
        setNu(.x) <- ifelse(1L == .y, value[[1]], value[[2]])
        .x
      })
    }

    invisible(object)
  })

setMethod("getRho", "H2ExtMO3FParam",
  function(object) {
    alpha <- getAlpha(object)

    3 * alpha / (4 - alpha)
  })

setMethod("getTau", "H2ExtMO3FParam",
  function(object) {
    alpha <- getAlpha(object)

    alpha / (2 - alpha)
  })

#' @importFrom purrr map_dbl
setMethod("getAlpha", "H2ExtMO3FParam",
  function(object) {
    fraction <- getFraction(object)
    alpha0 <- map_dbl(object@models, getAlpha)
    if (1L == length(alpha0)) {
      alpha <- alpha0 * fraction
    } else {
      alpha <- cumsum(alpha0[1:2] * c(fraction, 1 - fraction))
    }

    alpha
  })


#' @importFrom checkmate qassert
setValidity("H2ExtMO3FParam",
  function(object) {
    qassert(object@lambda, "N1(0,)")
    qassert(object@nu, "N2(0,)")

    invisible(TRUE)
  })


#' @importFrom purrr imap
setMethod("initialize", "H2ExtMO3FParam", # nolint
  function(.Object, # nolint
    partition = list(1L:2L, 3L:5L), lambda = 1e-1, nu = c(0.2, 0.3),
    fraction = 0.5, rho = NULL, tau = NULL, alpha = NULL) {
  if (!missing(partition) && !missing(lambda) &&
        (!missing(nu) || !missing(rho) || !missing(tau) || !missing(alpha)) &&
        !missing(fraction)) {
    .Object@fraction <- fraction
    if (missing(nu)) {
      if (!is.null(rho)) {
        nu <- invRho(.Object, rho)
      } else if (!is.null(tau)) {
        nu <- invTau(.Object, tau)
      } else if (!is.null(alpha)) {
        nu <- invAlpha(.Object, alpha)
      }
    }

    dim <- length(unlist(partition))
    models <- imap(c(dim, map_int(partition, length)), ~{
      new(getModelName(.Object),
            dim = .x, lambda = lambda, nu = nu[[pmin(.y, 2L)]])
      })

    .Object@dim <- dim
    .Object@partition <- partition
    .Object@models <- models
    .Object@lambda <- lambda
    .Object@nu <- nu

    validObject(.Object)
  }

  invisible(.Object)
  })


#' @rdname H2ExtMO3FParam-class
#'
#' @section Cuadras-Augé calibration parameter class:
#' Corresponds to a Lévy subordinators which are a convex combination of
#' a pure-killing subordinator and a pure-drift subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = \nu + (1 - \nu) x}
#'   \item \eqn{\alpha = \nu}
#' }
#'
#' @export CuadrasAugeH2ExtMO3FParam
CuadrasAugeH2ExtMO3FParam <- setClass("CuadrasAugeH2ExtMO3FParam", # nolint
  contains = "H2ExtMO3FParam")


#' @rdname H2ExtMO3FParam-class
#'
#' @section Alpha-stable calibration parameter class:
#' Corresponds to \eqn{\alpha}-stable subordinators.
#' \itemize{
#'   \item \eqn{\psi(x) = x^\nu}
#'   \item \eqn{\nu = \log_2(2 - \alpha)} and \eqn{\alpha = 2 - 2^\nu}
#' }
#'
#' @export AlphaStableH2ExtMO3FParam
AlphaStableH2ExtMO3FParam <- setClass("AlphaStableH2ExtMO3FParam", # nolint
  contains = "H2ExtMO3FParam")


#' @rdname H2ExtMO3FParam-class
#'
#' @section Poisson calibration parameter class:
#' Corresponds to Lévy subordinators which are a convex combination of
#' Poisson subordinators with jump size `nu` and pure-drift subordinators.
#' \itemize{
#'   \item \eqn{\psi(x) = \operatorname{e}^{-\nu}x + (1 - \operatorname{e}^{-x \nu})}
#'  \item \eqn{\nu = -log(1 - sqrt(\alpha))} and \eqn{\alpha = (1 - \operatorname{e}^{-\eta})}
#' }
#'
#' @export PoissonH2ExtMO3FParam
PoissonH2ExtMO3FParam <- setClass("PoissonH2ExtMO3FParam", # nolint
  contains = "H2ExtMO3FParam")


#' @rdname H2ExtMO3FParam-class
#'
#' @section Exponential calibration parameter class:
#' Corresponds to Lévy subordinator which is a convex combination of
#' Exponential-jump compound Poisson processes with rate `nu` and unit-intensity
#' and a pure-drift subordinators.
#' \itemize{
#'   \item \eqn{\psi(x) = (1 - 1 / (1 + \nu))x + 1 / (x + \nu)}
#'   \item \eqn{\nu = 0.5 \cdot (-3 + \sqrt{1 + 8 / \alpha})}
#'     and \eqn{\alpha = 2 / (1 + \nu) - 1 / (2 + \nu)}
#' }
#'
#' @export ExponentialH2ExtMO3FParam
ExponentialH2ExtMO3FParam <- setClass("ExponentialH2ExtMO3FParam", # nolint
  contains = "H2ExtMO3FParam")
