#' @include s4-H2ExCalibrationParam.R s4-H2ExtMO3FParam.R checkmate.R
NULL

#' Three-factor H2-extendible Gaussian calibration parameter
#'
#' [CalibrationParam-class] for the H2-extendible Gaussian copula with Exponential margin model for
#' the *(average) default counting process*  with 3 parameter. Extends [H2ExCalibration-class] and
#' related to [ExtGaussian2FParam-class].
#'
#' @slot lambda A non-negative number for the marginal rate.
#' @slot nu A numeric vector of length 2 for the model specific dependence parameters (global and
#'   component specific; range depends on specific model). Use `rho` or `tau` in the constructor to
#'   set dependence parameter.
#'
#' @details
#' The model is defined by the assumptoin that the *multivariate default times* \eqn{\tau = (\tau_1,
#' \ldots, \tau_d)} are from a H2-extendible Gaussian copula model with Exponential margins.
#' The model is specified by three parameters (in addition to the composition): The *marginal rate*
#' `lambda` and the (internal) *outer* and *inner dependency parameters* `nu` (Pearson correlation).
#' The dependence parameter `nu` should not be set by the user; instead they should provide either
#' `rho` (*Spearman's Rho*) or `tau` (*Kendall's Tau*).
#' The parameters `rho` or `tau` should be between zero and one, of length 2, and non-decreasing;
#' the first value represents the *outer dependence* between components of different partition
#' elements and the second value represents the *inner depenence* between components of the same
#' partition element.
#' The link between Spearman's Rho or Kendall's Tau and the
#' internal dependence parameter (Pearson correlation) is
#' \itemize{
#'   \item \eqn{\rho = 2 \sin(\rho_S \cdot \pi / 6)} and
#'     \eqn{\rho_S = 6 / \pi \cdot \arcsin(\rho/2)}
#'   \item \eqn{\rho = \sin(\tau \cdot \pi / 2)} and
#'     \eqn{\tau = 2 / \pi \cdot \arcsin(\rho)}
#' }
#'
#' For details on the underlying extendible model, see [ExtGaussian2FParam-class].
#'
#' @export H2ExtGaussian3FParam
H2ExtGaussian3FParam <- setClass("H2ExtGaussian3FParam", # nolint
  contains = "H2ExCalibrationParam",
  slots = c(lambda = "numeric", nu = "numeric"))


setMethod("getModelName", "H2ExtGaussian3FParam",
  function(object) {
    "ExtGaussian2FParam"
  })

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

setGeneric("calcNu2Rho",
  function(object, value) {
    standardGeneric("calcNu2Rho")
  })
#' @importFrom checkmate assert_numeric
setMethod("calcNu2Rho", "H2ExtGaussian3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    (6 / pi) * asin(value / 2)
  })

setMethod("getRho", "H2ExtGaussian3FParam",
  function(object) {
    getNu(object) %>%
      calcNu2Rho(object, .)
  })

#' @importFrom checkmate assert_numeric
setMethod("calcRho2Nu", "H2ExtGaussian3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    2 * sin((pi / 6) * value)
  })

#' @importFrom checkmate assert_numeric
setMethod("invRho", "H2ExtGaussian3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    calcRho2Nu(object, value)
  })

#' @importFrom checkmate assert_numeric
setReplaceMethod("setRho", "H2ExtGaussian3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    setNu(object) <- invRho(object, value)

    invisible(object)
  })

setGeneric("calcNu2Tau",
  function(object, value) {
    standardGeneric("calcNu2Tau")
  })
#' @importFrom checkmate assert_numeric
setMethod("calcNu2Tau", "H2ExtGaussian3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    (2 / pi) * asin(value)
  })

setMethod("getTau", "H2ExtGaussian3FParam",
  function(object) {
    getNu(object) %>%
      calcNu2Tau(object, .)
  })

#' @importFrom checkmate assert_numeric
setMethod("calcTau2Nu", "H2ExtGaussian3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    sin((pi / 2) * value)
  })

#' @importFrom checkmate assert_numeric
setMethod("invTau", "H2ExtGaussian3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    calcTau2Nu(object, value)
  })

setReplaceMethod("setTau", "H2ExtGaussian3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    setNu(object) <- invTau(object, value)

    invisible(object)
  })


#' @importFrom checkmate qassert assert_numeric
setValidity("H2ExtGaussian3FParam",
  function(object) {
    qassert(object@lambda, "N1(0,)")
    assert_numeric(object@nu, lower = 0, upper = 1, any.missing = TRUE, len = 2L, sorted = TRUE)

    invisible(TRUE)
  })

#' @describeIn H2ExtGaussian3FParam-class Constructor
#' @aliases initialize,H2ExtGaussian3FParam-method
#' @aliases initialize,H2ExtGaussian3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param composition An integerish vector with the composition.
#' @param lambda Marginal intensity.
#' @param rho *Outer* and *inner* bivariate Spearman's Rho.
#' @param tau *Outer* and *inner* bivariate Kendall's Tau.
#' @param nu (Internal) *Outer* and *inner* bivariate dependence parameter.
#'
#' @examples
#' H2ExtGaussian3FParam(composition = c(2L, 4L, 2L), lambda = 8e-2, rho = c(3e-1, 5e-1))
#' H2ExtGaussian3FParam(composition = c(2L, 4L, 2L), lambda = 8e-2, tau = c(3e-1, 5e-1))
setMethod("initialize", "H2ExtGaussian3FParam",
  function(.Object, # nolint
      composition = c(2L, 3L), lambda = 1e-1, nu = c(0.2, 0.3),
      rho = NULL, tau = NULL) {
    if (!missing(composition) && !missing(lambda) &&
          (!missing(nu) || !missing(rho) || !missing(tau))) {
      setComposition(.Object) <- composition
      setLambda(.Object) <- lambda
      if (!missing(nu)) {
        setNu(.Object) <- nu
      } else if (!missing(rho)) {
        setRho(.Object) <- rho
      } else if (!missing(tau)) {
        setTau(.Object) <- tau
      }

      validObject(.Object)
    }

    invisible(.Object)
  })


#' @describeIn H2ExtGaussian3FParam-class
#'    simulates the vector of *default times* and returns a matrix `x` with
#'    `dim(x) == c(n_sim, getDimension(object))`.
#' @aliases simulate_dt,H2ExtGaussian3FParam-method
#'
#' @inheritParams simulate_dt
#' @param n_sim Number of samples.
#'
#' @section Simulation:
#' The default times are sampled in a two-stage procedure: First a sample is drawn from the Gaussian
#' copula whose correlation matrix reflect the inner- and outer-dependency parameters, i.e.
#' \eqn{\rho_{i j} = \nu_1} if \eqn{i,j} are from different elements of the partition and
#' \eqn{\rho_{i j} = \nu_2} if \eqn{i j} are from the same element of the partition; then the
#' results are transformed using [stats::qexp()].
#'
#' @examples
#' parm <- H2ExtGaussian3FParam(composition = c(2L, 4L, 2L), lambda = 8e-2, rho = c(2e-1, 7e-1))
#' simulate_dt(parm, n_sim = 5L)
#'
#' @importFrom copula normalCopula rCopula P2p p2P
#' @importFrom stats qexp
#' @include utils.R
setMethod("simulate_dt", "H2ExtGaussian3FParam",
  function(object, ..., n_sim = 1e4) {
    d <- getDimension(object)
    lambda <- getLambda(object)
    nu <- getNu(object)
    corr <- p2P(nu[[1]], d = d)
    for (elem in getPartition(object)) {
      corr[elem, elem] <- p2P(nu[[2]], d = length(elem))
    }
    cop <- normalCopula(param = P2p(corr), dim = d, dispstr = "un")

    qexp(rCopula(n_sim, cop), rate = lambda, lower.tail = FALSE)
  })


#' @describeIn H2ExtGaussian3FParam-class Display the object.
#' @aliases show,H2ExtGaussian3FParam-method
#'
#' @inheritParams methods::show
#'
#' @export
setMethod("show", "H2ExtGaussian3FParam",
 function(object) {
   cat(sprintf("An object of class %s\n", classLabel(class(object))))
   cat(sprintf("Partition: %s = %s\n", getDimension(object),
     paste(getComposition(object), collapse = " + ")))
   to_vector <- function(x) {
     paste0("(", paste(x, collapse = ", "), ")")
   }
   cat("Parameter:\n")
   cat(sprintf("* %s: %s\n", "Lambda", format(getLambda(object))))
   cat(sprintf("* %s: %s\n", "Rho", to_vector(format(getRho(object)))))
   cat(sprintf("* %s: %s\n", "Tau", to_vector(format(getTau(object)))))
   cat("Internal parameter:\n")
   cat(sprintf("* %s: %s\n", "Nu", to_vector(format(getNu(object)))))
  })
