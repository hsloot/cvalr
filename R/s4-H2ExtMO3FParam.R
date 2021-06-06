#' @include s4-H2ExtMOParam.R s4-ExtMO2FParam.R checkmate.R
NULL

# nolint start
ERR_MSG_NU2 <- "`nu` must be scalar double"
ERR_MSG_NU2_INTERVAL <- paste(ERR_MSG_NU2, "in interval %s^2")
# nolint end

#' Three-factor H2-extendible Marshall--Olkin calibration parameter
#'
#' [CalibrationParam-class] for the H2-extendible Marshall-Olkin *(average) default counting
#' process* model with 3 parameter. Extends [H2ExtMOParam-class] and related to
#' [ExtMO2FParam-class].
#'
#' @slot lambda A non-negative number for the marginal rate.
#' @slot nu A numeric vector of length 2 for the model specific dependence parameters (global and
#'   component specific; range depends on specific model). Use `rho`, `tau`, or `alpha` in the
#'   constructor to set dependence parameter.
#'
#' @details
#' The model is defined by the assumption that the *multivariate default times* \eqn{\tau = (\tau_1,
#' \ldots, \tau_d)} are H2-extendible Marshall-Olkin, see [H2ExtMOParam-class] for the details. This
#' class provides an interface for easy-to-use, 3-factor families for this model. For all
#' implemented families, the *marginal rate* can be specified by `lambda` and the *(internal)
#' dependence parameters* (model specific) of the global model and the component models can be
#' specified by `nu`. The dependence parameter `nu` should not be set by the user; instead they
#' should provide either `rho` (*Spearman's Rho*), `tau` (*Kendall's Tau*), or `alpha` (*lower
#' tail-dependence coefficient*).
#' The parameters `rho`, `tau`, or `alpha` should be between zero and one, of length 2, and
#' non-decreasing; the first value represents the *outer dependence* between components of different
#' partition elements and the second value represents the *inner dependence* between components of
#' the same partition element. Setting either of the three dependence parameters implicitely sets
#' the `fraction`-slot, too.
#' The link between lower tail-dependence coefficient \eqn{\alpha} and
#' Spearman's Rho and Kendall's Tau is (all calculations are component-wise)
#' \itemize{
#'   \item \eqn{\alpha = 4 \rho / (3 + \rho)} and \eqn{\rho = 3 \alpha / (4 - \alpha)}
#'   \item \eqn{\alpha = 2 \tau / (1 + \tau)} and \eqn{\tau = \alpha / (2 - \alpha)}
#' }
#' Consider \eqn{\tilde{\alpha}} to be the vector of the actual lower TDC of the global and the
#' component model and \eqn{\kappa} to be the `fraction` parameter. Then:
#' \deqn{
#'   \alpha_1 = \kappa \tilde{\alpha}_1,
#'     \alpha_2 = \kappa \tilde{\alpha}_1 + (1 - \kappa) \tilde{\alpha}_2
#' }
#' and
#' \deqn{
#'   \tilde{\alpha}_1 = \frac{\alpha_1}{\kappa},
#'     \tilde{\alpha_2} = \frac{\alpha_2 - \alpha_1}{1 - \kappa} .
#' }
#' In particular, the boundaries \eqn{\tilde{\alpha}_i \in [0, 1]} impose the restrictions
#' \deqn{
#'   \alpha_1 \leq \kappa \leq \alpha_1 + (1 - \alpha_2).
#' }
#' For the families deriving from [H2ExtMO3FParam-class] we choose the default value for fraction to
#' be the midpoint of the admissible interval, i.e.
#' \deqn{
#'   \kappa
#'     = \frac{2 \alpha_1 + 1 - \alpha_2}{2} .
#' }
#'
#' For details on the underlying extendible models, see [ExtMO2FParam-class].
#'
#' @export
setClass("H2ExtMO3FParam", # nolint
  contains = c("H2ExtMOParam", "VIRTUAL"),
  slots = c(lambda = "numeric", nu = "numeric"))


setGeneric("constructModels",
  function(object, ...) {
    standardGeneric("constructModels")
  })
#' @importFrom checkmate test_list qtest test_integerish test_numeric
#' @importFrom purrr map imap
setMethod("constructModels", "H2ExtMO3FParam",
function(object, composition = getComposition(object), lambda = getLambda(object),
         nu = getNu(object)) {
  dim <- sum(composition)
  if (test_integerish(composition, lower = 1L, any.missing = FALSE, min.len = 1L) &&
      qtest(lambda, "N1(0,)") && test_numeric(nu, any.missing = FALSE, len = 2L)) {
    out <- imap(c(dim, composition), ~{
      new(getModelName(object), .x, lambda, nu[[pmin(.y, 2L)]])
      })
  } else {
    out <- list()
  }

  out
})

#' @importFrom rmo exIntensities
#' @importFrom checkmate test_list qassert
setReplaceMethod("setComposition", "H2ExtMO3FParam",
  function(object, value) {
    qassert(value, "X+[1,)")
    dim <- sum(value)
    qassert(dim, "X1[2,)")
    if (isTRUE(all(getComposition(object) != value))) {
      object <- callNextMethod()
      models <- constructModels(object)
      if (test_list(models, types = getModelName(object), any.missing = FALSE, min.len = 2L)) {
        setModels(object) <- models
      }
    }

    invisible(object)
  })

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
    object@lambda <- value
    models <- constructModels(object)
    if (test_list(models, types = getModelName(object), any.missing = FALSE, min.len = 2L)) {
      setModels(object) <- models
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
    object@nu <- value
    models <- constructModels(object)
    if (test_list(models, types = getModelName(object), any.missing = FALSE, min.len = 2L)) {
      setModels(object) <- models
    }

    invisible(object)
  })

#' @importFrom checkmate qassert
setMethod("calcAlpha2Rho", "H2ExtMO3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    3 * value / (4 - value)
  })

setMethod("getRho", "H2ExtMO3FParam",
  function(object) {
    getAlpha(object) %>%
      calcAlpha2Rho(object, .)
  })

#' @importFrom checkmate assert_numeric
setMethod("calcAlpha2Tau", "H2ExtMO3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    value / (2 - value)
  })

setMethod("getTau", "H2ExtMO3FParam",
  function(object) {
    getAlpha(object) %>%
      calcAlpha2Tau(object, .)
  })

#' @importFrom utils head
#' @importFrom purrr map_dbl
setMethod("getAlpha", "H2ExtMO3FParam",
  function(object) {
    fraction <- getFraction(object)
    map_dbl(getModels(object), getAlpha) %>%
      head(2L) %>%
      `*`(., c(fraction, 1 - fraction)) %>%
      cumsum
  })

#' @importFrom checkmate assert_numeric
setMethod("calcRho2Alpha", "H2ExtMO3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    4 * value / (3 + value)
  })

#' @importFrom checkmate assert_numeric
setMethod("calcTau2Alpha", "H2ExtMO3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    2 * value / (1 + value)
  })

setGeneric("calcAlpha2Fraction",
  function(object, value) {
    standardGeneric("calcAlpha2Fraction")
  })
#' @importFrom checkmate test_numeric
setMethod("calcAlpha2Fraction", "H2ExtMO3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    value[[1]] + 0.5  * (1 - value[[2]])
  })

#' @importFrom checkmate test_numeric
setReplaceMethod("setAlpha", "H2ExtMO3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    setFraction(object) <- calcAlpha2Fraction(object, value)
    setNu(object) <- invAlpha(object, value)

    invisible(object)
  })

#' @importFrom checkmate test_numeric
setReplaceMethod("setRho", "H2ExtMO3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    setAlpha(object) <- calcRho2Alpha(object, value)

    invisible(object)
  })

#' @importFrom checkmate test_numeric
setReplaceMethod("setTau", "H2ExtMO3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    setAlpha(object) <- calcTau2Alpha(object, value)

    invisible(object)
  })



#' @importFrom checkmate qtest test_numeric
setValidity("H2ExtMO3FParam",
  function(object) {
    if (!qtest(object@lambda, "N1(0,)")) {
      return(ERR_MSG_LAMBDA)
    }
    if (!qtest(object@nu, "N2(0,)")) {
      return(sprintf(ERR_MSG_NU2_INTERVAL, "(0,)"))
    }

    invisible(TRUE)
  })

#' @describeIn H2ExtMO3FParam-class Constructor
#' @aliases initialize,H2ExtMO3FParam-method
#' @aliases initialize,H2ExtMO3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param composition An integerish vector with the composition.
#' @param lambda Marginal intensity.
#' @param rho *Outer* and *inner* bivariate Spearman's Rho.
#' @param tau *Outer* and *inner* bivariate Kendall's Tau.
#' @param alpha *Outer* and *inner* bivariate lower tail-dependence coefficient.
#' @param nu (Internal) *Outer* and *inner* bivariate dependence parameter.
#' @param fraction (Internal) proportion associated with the global model, see details.
#'
#' @examples
#' CuadrasAugeH2ExtMO3FParam()
#' CuadrasAugeH2ExtMO3FParam(composition = c(2L, 4L, 2L), lambda = 8e-2, rho = c(3e-1, 5e-1))
#' AlphaStableH2ExtMO3FParam()
#' AlphaStableH2ExtMO3FParam(composition = c(2L, 4L, 2L), lambda = 8e-2, rho = c(3e-1, 5e-1))
#' PoissonH2ExtMO3FParam()
#' PoissonH2ExtMO3FParam(composition = c(2L, 4L, 2L), lambda = 8e-2, tau = c(3e-1, 5e-1))
#' ExponentialH2ExtMO3FParam()
#' ExponentialH2ExtMO3FParam(composition = c(2L, 4L, 2L), lambda = 8e-2, alpha = c(3e-1, 5e-1))
setMethod("initialize", "H2ExtMO3FParam", # nolint
  function(.Object, # nolint
    composition, lambda, nu, fraction, rho, tau, alpha) {
  if (!missing(composition) && !missing(lambda) &&
        ((!missing(nu) && !missing(fraction)) || !missing(rho) || !missing(tau) ||
         !missing(alpha))) {
    setComposition(.Object) <- composition
    setLambda(.Object) <- lambda
    if (!missing(nu)) {
      setFraction(.Object) <- ifelse(missing(fraction), 0.5, fraction)
      setNu(.Object) <- nu
    } else if (!missing(rho)) {
      setRho(.Object) <- rho
    } else if (!missing(tau)) {
      setTau(.Object) <- tau
    } else if (!missing(alpha)) {
      setAlpha(.Object) <- alpha
    }

    validObject(.Object)
  }

  invisible(.Object)
  })


setMethod("getModelName", "H2ExtMO3FParam",
  function(object) {
    "ExtMO2FParam"
  })

#' @describeIn H2ExtMO3FParam-class
#'   calculates the *payoff equation* for a *portfolio CDS* (vectorized w.r.t.
#'   the argumentes `recovery_rate`, `coupon`, and `upfront`).
#' @aliases expected_pcds_equation,H2ExtMO3FParam-method
#'
#' @inheritParams expected_pcds_equation
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#'
#' @inheritSection ExtMO2FParam-class Expected portfolio CDS loss
#'
#' @examples
#' parm <- AlphaStableH2ExtMO3FParam(c(3, 3, 4, 5), 8e-2, rho = c(3e-1, 6e-1))
#' expected_pcds_equation(
#'   parm, times = seq(0.25, 5, by = 0.25), discount_factors = rep(1, 20L), recovery_rate = 0.4,
#'   coupon = 1e-1, upfront = 0)
#' expected_pcds_equation(
#'   parm, times = seq(0.25, 5, by = 0.25), discount_factors = rep(1, 20L), recovery_rate = 0.4,
#'   coupon = 1e-1, upfront = 0, method = "mc", n_sim = 1e1)
#'
#' @importFrom stats pexp
#' @importFrom checkmate assert check_numeric
#' @importFrom vctrs vec_size_common vec_recycle
#' @include RcppExports.R
#'
#' @export
setMethod("expected_pcds_equation", "H2ExtMO3FParam",
  function(object, times, discount_factors, recovery_rate, coupon, upfront, ...,
      method = c("default", "prob", "mc")) {
    method <- match.arg(method)
    if ("default" == method) {
      qassert(times, "N+[0,)")
      qassert(discount_factors, paste0("N", length(times), "[0,)"))
      qassert(recovery_rate, "N+[0,1]")
      qassert(coupon, "N+")
      qassert(upfront, "N+")

      p <- vec_size_common(recovery_rate, coupon, upfront)
      recovery_rate <- vec_recycle(recovery_rate, p)
      coupon <- vec_recycle(coupon, p)
      upfront <- vec_recycle(upfront, p)

      x <- pexp(times, rate = getLambda(object)) %*% t(1 - recovery_rate)

      out <- Rcpp__lagg_ev_pcds(x, times, discount_factors, recovery_rate, coupon, upfront)
    } else {
      out <- callNextMethod(
        object, times, discount_factors, recovery_rate, coupon, upfront, ..., method = method)
    }

    out
  })


#' @describeIn H2ExtMO3FParam-class Display the object.
#' @aliases show,H2ExtMO3FParam-method
#'
#' @param object A [CalibrationParam-class]-object.
#'
#' @importFrom utils capture.output
#' @importFrom purrr map compose flatten_chr
#'
#' @export
setMethod("show", "H2ExtMO3FParam",
  function(object) {
    cat(sprintf("An object of class %s\n", classLabel(class(object))))
    if (isTRUE(validObject(object, test = TRUE))) {
      cat(sprintf("Composition: %s = %s\n", getDimension(object),
        paste(getComposition(object), collapse = " + ")))
      to_vector <- function(x) {
        paste0("(", paste(x, collapse = ", "), ")")
      }
      cat("Parameter:\n")
      cat(sprintf("* %s: %s\n", "Lambda", format(getLambda(object))))
      cat(sprintf("* %s: %s\n", "Rho", to_vector(format(getRho(object)))))
      cat(sprintf("* %s: %s\n", "Tau", to_vector(format(getTau(object)))))
      cat(sprintf("* %s: %s\n", "Alpha", to_vector(format(getAlpha(object)))))
      cat("Internal parameter:\n")
      cat(sprintf("* %s: %s\n", "Nu", to_vector(format(getNu(object)))))
      cat(sprintf("* %s: %s\n", "Fraction", format(getFraction(object))))
      cat("Models:\n")
      cat("* Global model\n")
      writeLines(
        paste0("\t", capture.output(show(as(getGlobalModel(object), getModelName(object))))))
      cat("* Partition models:\n")
      to_list_item <- function(x) {
        out <- rep("  ", length(x))
        out[[1]] <- "- "

        paste0(out, x)
      }
      getPartitionModels(object) %>%
        map(compose(to_list_item, ~capture.output(show(.)), ~as(., getModelName(object)))) %>%
        flatten_chr %>%
        paste0("\t", .) %>%
        writeLines
    } else {
      cat("\t (invalid or not initialized)\n")
    }

    invisible(NULL)
    })



#' @rdname H2ExtMO3FParam-class
#'
#' @inheritSection ExtMO2FParam-class Cuadras-Augé calibration parameter class
#'
#' @export CuadrasAugeH2ExtMO3FParam
CuadrasAugeH2ExtMO3FParam <- setClass("CuadrasAugeH2ExtMO3FParam", # nolint
  contains = "H2ExtMO3FParam")


setMethod("getModelName", "CuadrasAugeH2ExtMO3FParam",
  function(object) {
    "CuadrasAugeExtMO2FParam"
  })

#' @include utils.R
#' @importFrom checkmate qassert
setMethod("invAlpha", "CuadrasAugeH2ExtMO3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    fraction <- getFraction(object)
    adjacent_differences(value) / c(fraction, 1 - fraction)
  })



#' @rdname H2ExtMO3FParam-class
#'
#' @inheritSection ExtMO2FParam-class Alpha-stable calibration parameter class
#'
#' @export AlphaStableH2ExtMO3FParam
AlphaStableH2ExtMO3FParam <- setClass("AlphaStableH2ExtMO3FParam", # nolint
  contains = "H2ExtMO3FParam")


setMethod("getModelName", "AlphaStableH2ExtMO3FParam",
  function(object) {
    "AlphaStableExtMO2FParam"
  })

#' @include utils.R
#' @importFrom checkmate qassert
setMethod("invAlpha", "AlphaStableH2ExtMO3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    fraction <- getFraction(object)
    value <- adjacent_differences(value) / c(fraction, 1 - fraction)
    log2(2 - value)

  })



#' @rdname H2ExtMO3FParam-class
#'
#' @inheritSection ExtMO2FParam-class Poisson calibration parameter class
#'
#' @export PoissonH2ExtMO3FParam
PoissonH2ExtMO3FParam <- setClass("PoissonH2ExtMO3FParam", # nolint
  contains = "H2ExtMO3FParam")


setMethod("getModelName", "PoissonH2ExtMO3FParam",
  function(object) {
    "PoissonExtMO2FParam"
  })

#' @include utils.R
#' @importFrom checkmate qassert
setMethod("invAlpha", "PoissonH2ExtMO3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    fraction <- getFraction(object)
    value <- adjacent_differences(value) / c(fraction, 1 - fraction)
    -log(1 - sqrt(value))
  })



#' @rdname H2ExtMO3FParam-class
#'
#' @inheritSection ExtMO2FParam-class Exponential calibration parameter class
#'
#' @export ExponentialH2ExtMO3FParam
ExponentialH2ExtMO3FParam <- setClass("ExponentialH2ExtMO3FParam", # nolint
  contains = "H2ExtMO3FParam")


setMethod("getModelName", "ExponentialH2ExtMO3FParam",
  function(object) {
    "ExponentialExtMO2FParam"
  })

#' @include utils.R
#' @importFrom checkmate qassert
setMethod("invAlpha", "ExponentialH2ExtMO3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    fraction <- getFraction(object)
    value <- adjacent_differences(value) / c(fraction, 1 - fraction)
    0.5 * (-3 + sqrt(1 + 8 / value))
  })
