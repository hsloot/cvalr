#' @include s4-ExtMOParam.R checkmate.R
NULL

# nolint start
ERR_MSG_LAMBDA <- "`lambda` must be positive scalar double"
ERR_MSG_NU1 <- "`nu` must be scalar double"
ERR_MSG_NU1_INTERVAL <- paste(ERR_MSG_NU1, "in interval %s")
ERR_MSG_BF_TYPE <- "`bf` is of wrong type"
# nolint end

#' Two-factor extendible Marshall-Olkin calibration parameters
#'
#' @description
#' [CalibrationParam-class]-class with two parameters for the extendible
#' Marshall-Olkin model for the *(average) default counting process*.
#' Extends [ExtMOParam-class].
#'
#' @slot lambda A non-negative number for the marginal rate.
#' @slot nu A numeric number for the model specific dependence parameter (range
#'   depends on specific model, use `rho`, `tau`, or `alpha` to set dependence
#'   parameter).
#'
#' @details
#' The model is defined by the assumption that the multivariate default times
#' \eqn{\tau = (\tau_1, \ldots, \tau_d)} are extendible Marshall-Olkin, see
#' [ExtMOParam-class] for the details. This class provides an interface for
#' easy-to-use, 2-factor families for this model.
#' For all implemented families, the marginal rate can be specified by `lambda`
#' and the dependence can be specified by the internal parameter `nu`.
#' For all implemented families, the (internal) dependence parameter `nu` has a
#' one-to-one relationship, and can be replaced by, *Spearman's Rho* `rho`,
#' *Kendall' Tau* `tau` or the *(lower) tail dependence coefficient* `alpha`.
#' The possible range for `rho`, `tau`, and `alpha` is from zero to one
#' (boundaries might not be included).
#' The link between lower tail-dependence coefficient \eqn{\alpha} and
#' Spearman's Rho and Kendall's Tau is
#' \itemize{
#'   \item \eqn{\alpha = 4 \rho / (3 + \rho)} and \eqn{\rho = 3 \alpha / (4 - \alpha)}
#'   \item \eqn{\alpha = 2 \tau / (1 + \tau)} and \eqn{\tau = \alpha / (2 - \alpha)}
#' }
#'
#' @export
setClass("ExtMO2FParam", # nolint
  contains = c("ExtMOParam", "VIRTUAL"),
  slots = c(lambda = "numeric", nu = "numeric"))


setGeneric("constructBernsteinFunction",
  function(object, ...) {
    standardGeneric("constructBernsteinFunction")
  })

#' @importFrom rmo ScaledBernsteinFunction valueOf
#' @importFrom checkmate assert check_class
setReplaceMethod("setBernsteinFunction", "ExtMO2FParam",
  function(object, value) {
    assert(combine = "and",
      check_class(value, "ScaledBernsteinFunction"),
      check_equal(1, valueOf(value@original, 1, 0L)))

    callNextMethod()
  })

setGeneric("getLambda",
  function(object) {
    standardGeneric("getLambda")
  })
setMethod("getLambda", "ExtMO2FParam",
  function(object) {
    object@lambda
  })

setGeneric("setLambda<-",
  function(object, value) {
    standardGeneric("setLambda<-")
  })
#' @importFrom checkmate qassert qtest
setReplaceMethod("setLambda", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    object@lambda <- value
    nu <- getNu(object)
    if (qtest(nu, "N1")) {
      setBernsteinFunction(object) <- constructBernsteinFunction(object, value, nu)
    }

    invisible(object)
  })

setGeneric("getNu",
  function(object) {
    standardGeneric("getNu")
  })
setMethod("getNu", "ExtMO2FParam",
  function(object) {
    object@nu
  })

setGeneric("setNu<-",
  function(object, value) {
    standardGeneric("setNu<-")
  })
#' @importFrom checkmate qassert qtest
setReplaceMethod("setNu", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1")
    object@nu <- value
    lambda <- getLambda(object)
    if (qtest(lambda, "N1(0,)")) {
      setBernsteinFunction(object) <- constructBernsteinFunction(object, lambda, value)
    }

    invisible(object)
  })

setGeneric("calcAlpha2Nu",
  function(object, value) {
    standardGeneric("calcAlpha2Nu")
  })

setGeneric("calcNu2Alpha",
  function(object, value) {
    standardGeneric("calcNu2Alpha")
  })

#' @importFrom checkmate qassert
setGeneric("invAlpha",
  function(object, value) {
    standardGeneric("invAlpha")
  })
#' @importFrom checkmate qassert
setMethod("invAlpha", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    calcAlpha2Nu(object, value)
  })

setGeneric("getAlpha",
  function(object) {
    standardGeneric("getAlpha")
  })
#' @importFrom rmo valueOf
setMethod("getAlpha", "ExtMO2FParam",
  function(object) {
    getNu(object) %>%
      calcNu2Alpha(object, .)
  })

setGeneric("setAlpha<-",
  function(object, value) {
    standardGeneric("setAlpha<-")
  })
#' @importFrom checkmate qassert
setReplaceMethod("setAlpha", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invAlpha(object, value)

    invisible(object)
  })

setGeneric("calcRho2Alpha",
  function(object, value) {
    standardGeneric("calcRho2Alpha")
  })
#' @importFrom checkmate qassert
setMethod("calcRho2Alpha", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    4 * value / (3 + value)
  })

setGeneric("calcRho2Nu",
  function(object, value) {
    standardGeneric("calcRho2Nu")
  })
#' @importFrom checkmate qassert
setMethod("calcRho2Nu", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    calcRho2Alpha(object, value) %>%
      calcAlpha2Nu(object, .)
  })

setGeneric("invRho",
  function(object, value) {
    standardGeneric("invRho")
  })
#' @importFrom checkmate qassert
setMethod("invRho", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    calcRho2Nu(object, value)
  })

setGeneric("calcTau2Alpha",
  function(object, value) {
    standardGeneric("calcTau2Alpha")
  })
#' @importFrom checkmate qassert
setMethod("calcTau2Alpha", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    2 * value / (1 + value)
  })

setGeneric("calcTau2Nu",
  function(object, value) {
    standardGeneric("calcTau2Nu")
  })
#' @importFrom checkmate qassert
setMethod("calcTau2Nu", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    calcTau2Alpha(object, value) %>%
      calcAlpha2Nu(object, .)
  })

setGeneric("invTau",
  function(object, value) {
    standardGeneric("invTau")
  })
#' @importFrom checkmate qassert
setMethod("invTau", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    calcTau2Nu(object, value)
  })

setGeneric("calcAlpha2Rho",
  function(object, value) {
    standardGeneric("calcAlpha2Rho")
  })
#' @importFrom checkmate qassert
setMethod("calcAlpha2Rho", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    3 * value / (4 - value)
  })

setGeneric("getRho",
  function(object) {
    standardGeneric("getRho")
  })
#' @importFrom checkmate qassert
setMethod("getRho", "ExtMO2FParam",
  function(object) {
    getAlpha(object) %>%
      calcAlpha2Rho(object, .)
  })

setGeneric("setRho<-",
  function(object, value) {
    standardGeneric("setRho<-")
  })
#' @importFrom checkmate qassert
setReplaceMethod("setRho", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invRho(object, value)

    invisible(object)
  })

setGeneric("calcAlpha2Tau",
  function(object, value) {
    standardGeneric("calcAlpha2Tau")
  })
#' @importFrom checkmate qassert
setMethod("calcAlpha2Tau", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    value / (2 - value)
  })

setGeneric("getTau",
  function(object) {
    standardGeneric("getTau")
  })
setMethod("getTau", "ExtMO2FParam",
  function(object) {
    getAlpha(object) %>%
      calcAlpha2Tau(object, .)
  })

setGeneric("setTau<-",
  function(object, value) {
    standardGeneric("setTau<-")
  })
#' @importFrom checkmate qassert
setReplaceMethod("setTau", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invTau(object, value)

    invisible(object)
  })



#' @importFrom rmo valueOf ScaledBernsteinFunction
#' @importFrom checkmate qtest test_class test_choice
#' @include checkmate.R
setValidity("ExtMO2FParam", # nolint
  function(object) {
    if (!qtest(object@lambda, "N1(0,)")) {
      return(ERR_MSG_LAMBDA)
    }
    if (!qtest(object@nu, "N1")) {
      return(ERR_MSG_NU1)
    }
    bf <- getBernsteinFunction(object)
    if (!(
        test_class(bf, "ScaledBernsteinFunction") && isTRUE(validObject(bf, test = TRUE)) &&
          test_choice(bf@scale, getLambda(object)) &&
          test_equal(1, valueOf(bf@original, 1, 0L)) &&
          test_equal(getNu(object), invAlpha(object, 2 - valueOf(bf@original, 2, 0L))))) {
      return(ERR_MSG_BF_TYPE)
    }

    invisible(TRUE)
  })


#' @describeIn ExtMO2FParam-class Constructor
#' @aliases initialize,ExtMO2FParam-method
#' @aliases initialize,ExtMO2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param dim Dimension.
#' @param lambda Marginal intensity.
#' @param nu (Internal) bivariate dependence parameter.
#' @param rho Bivariate Spearman's Rho.
#' @param tau Bivariate Kendall's Tau.
#' @param alpha Bivariate lower tail-dependence coefficient.
#'
#' @examples
#' ArmageddonExtMO2FParam()
#' ArmageddonExtMO2FParam(dim = 5L, lambda = 8e-2, rho = 4e-1)
#' AlphaStableExtMO2FParam()
#' AlphaStableExtMO2FParam(dim = 5L, lambda = 8e-2, rho = 4e-1)
#' PoissonExtMO2FParam()
#' PoissonExtMO2FParam(dim = 5L, lambda = 8e-2, tau = 4e-1)
#' ExponentialExtMO2FParam()
#' ExponentialExtMO2FParam(dim = 5L, lambda = 8e-2, alpha = 4e-1)
setMethod("initialize", signature = "ExtMO2FParam", # nolint
  definition = function(.Object, # nolint
      dim, lambda, nu, rho, tau, alpha) {
    if (!missing(dim) && !missing(lambda) &&
          !(missing(nu) && missing(rho) && missing(tau) && missing(alpha))) {
      setDimension(.Object) <- dim
      setLambda(.Object) <- lambda
      if (!missing(nu)) {
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



#' @describeIn ExtMO2FParam-class
#'   calculates the *payoff equation* for a *portfolio CDS* (vectorized w.r.t.
#'   the argumentes `recovery_rate`, `coupon`, and `upfront`).
#' @aliases expected_pcds_equation,ExtMO2FParam-method
#'
#' @inheritParams expected_pcds_equation
#'
#' @section Expected portfolio CDS loss:
#' The *expected portfolio CDS loss* for *recovery rate* \eqn{R} is calculated
#' using that
#' \deqn{
#'   \mathbb{E}[g(L_t)]
#'     = (1 - R) \cdot F(t)
#' }
#' with \eqn{g(x) = (1 - R) \cdot x} and \eqn{F} being the Exponential
#' distribution function for rate \eqn{\lambda}.
#'
#' @examples
#' parm <- ArmageddonExtMO2FParam(dim = 75L, lambda = 0.05, rho = 0.4)
#' expected_pcds_equation(
#'   parm, times = seq(25e-2, 5, by = 25e-2), discount_factors = rep(1, 20L), recovery_rate = 0.4,
#'   coupon = 1e-1, upfront = 0)
#' expected_pcds_equation(
#'   parm, times = seq(25e-2, 5, by = 25e-2), discount_factors = rep(1, 20L), recovery_rate = 0.4,
#'   coupon = 1e-1, upfront = 0, method = "mc", n_sim = 1e1)
#'
#' @importFrom stats pexp
#' @importFrom checkmate assert check_numeric
#' @importFrom vctrs vec_size_common vec_recycle
#' @include RcppExports.R
#'
#' @export
setMethod("expected_pcds_equation", "ExtMO2FParam",
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


#' @describeIn ExtMO2FParam-class Display the object.
#' @aliases show,ExtMO2FParam-method
#'
#' @inheritParams methods::show
#'
#' @export
setMethod("show", "ExtMO2FParam",
 function(object) {
   cat(sprintf("An object of class %s\n", classLabel(class(object))))
   if (isTRUE(validObject(object, test = TRUE))) {
     cat(sprintf("Dimension: %i\n", getDimension(object)))
     cat("Parameter:\n")
     cat(sprintf("* %s: %s\n", "Lambda", format(getLambda(object))))
     cat(sprintf("* %s: %s\n", "Rho", format(getRho(object))))
     cat(sprintf("* %s: %s\n", "Tau", format(getTau(object))))
     cat(sprintf("* %s: %s\n", "Alpha", format(getAlpha(object))))
     cat("Internal parameter:\n")
     cat(sprintf("* %s: %s\n", "Nu", format(getNu(object))))
   } else {
     cat("\t (invalid or not initialized)\n")
   }

   invisible(TRUE)
  })



#' @rdname ExtMO2FParam-class
#'
#' @section Armageddon-shock calibration parameter class:
#' Corresponds to a Lévy subordinator which is a convex combination of
#' a pure-killing subordinator and a pure-drift subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = \nu + (1 - \nu) x}
#'   \item \eqn{\alpha = \nu}
#' }
#'
#' @export ArmageddonExtMO2FParam
ArmageddonExtMO2FParam <- setClass("ArmageddonExtMO2FParam", # nolint
  contains = "ExtMO2FParam")

#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions LinearBernsteinFunction
#'   ConstantBernsteinFunction
#' @importFrom checkmate qtest test_class
setValidity("ArmageddonExtMO2FParam",
  function(object) {
    if (!qtest(object@nu, "N1[0,1]")) {
      return(sprintf(ERR_MSG_NU1_INTERVAL, "[0,1]"))
    }
    bf <- getBernsteinFunction(object)
    if (!(
        test_class(bf@original, "SumOfBernsteinFunctions") &&
          test_class(bf@original@first, "LinearBernsteinFunction") &&
          test_class(bf@original@second, "ConstantBernsteinFunction"))) {
      return(ERR_MSG_BF_TYPE)
    }

    invisible(TRUE)
  })


#' @importFrom rmo SumOfBernsteinFunctions LinearBernsteinFunction ConstantBernsteinFunction
#' @importFrom checkmate qassert
setMethod("constructBernsteinFunction", "ArmageddonExtMO2FParam",
  function(object, lambda, nu, ...) {
    qassert(lambda, "N1(0,)")
    qassert(nu, "N1[0,1]")
    ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = 1 - nu),
        second = ConstantBernsteinFunction(constant = nu))
    )
  })

#' @importFrom checkmate qassert
setMethod("calcAlpha2Nu", "ArmageddonExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    value
  })

#' @importFrom checkmate qassert
setMethod("calcNu2Alpha", "ArmageddonExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    value
  })

#' @importFrom checkmate qassert
setReplaceMethod("setNu", "ArmageddonExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    callNextMethod()
  })

#' @describeIn ExtMO2FParam-class
#'    simulates the vector of *default times* and returns a matrix `x` with
#'    `dim(x) == c(n_sim, getDimension(object))`.
#' @aliases simulate_dt,ArmageddonExtMO2FParam-method
#'
#' @inheritParams simulate_dt
#' @param n_sim Number of samples.
#'
#' @section Simulation:
#' The default times are sampled using [rmo::rarmextmo_esm()].
#'
#'
#' @examples
#' parm <- ArmageddonExtMO2FParam(dim = 5L, lambda = 8e-2, rho = 4e-1)
#' simulate_dt(parm, n_sim = 5L)
#'
#' @importFrom rmo rarmextmo_esm
#' @include utils.R
#'
#' @export
setMethod("simulate_dt", "ArmageddonExtMO2FParam",
  function(object, ..., n_sim = 1e4L) {
    rarmextmo_esm(
      n_sim, getDimension(object),
      getLambda(object) * (1 - getNu(object)), getLambda(object) * getNu(object))
  })


#' @describeIn ExtMO2FParam-class
#'   simulates the *average default counting process* and returns a
#'   matrix `x` with `dim(x) == c(n_sim, length(times))`.
#' @aliases simulate_adcp,ArmageddonExtMO2FParam-method
#'
#' @inheritParams simulate_adcp
#' @param times A non-negative numeric vector of timepoints.
#' @param n_sim Number of samples.
#'
#' @section Simulation:
#' The default times are sampled using [rmo::rexmo_markovian()].
#'
#'
#' @examples
#' parm <- ArmageddonExtMO2FParam(dim = 5L, lambda = 8e-2, rho = 4e-1)
#' simulate_adcp(parm, times = seq(25e-2, 5, by = 25e-2), n_sim = 5L)
#'
#' @include RcppExports.R
#'
#' @export
setMethod("simulate_adcp", "ArmageddonExtMO2FParam",
  function(object, times, ..., n_sim = 1e4L) {
    Rcpp__rarmextmo_esm_adcp(
      n_sim, times, getDimension(object),
      getLambda(object) * (1 - getNu(object)), getLambda(object) * getNu(object))
  })


#' @rdname ExtMO2FParam-class
#'
#' @section Alpha-stable calibration parameter class:
#' Corresponds to an \eqn{\alpha}-stable subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = x^\nu}
#'   \item \eqn{\nu = \log_2(2 - \alpha)} and \eqn{\alpha = 2 - 2^\nu}
#' }
#'
#' @export AlphaStableExtMO2FParam
AlphaStableExtMO2FParam <- setClass("AlphaStableExtMO2FParam", # nolint
  contains = "ExtMO2FParam")

#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions AlphaStableBernsteinFunction
#' @importFrom checkmate qtest test_class
setValidity("AlphaStableExtMO2FParam",
  function(object) {
    if (!qtest(object@nu, "N1(0,1)")) {
      return(sprintf(ERR_MSG_NU1_INTERVAL, classLabel(class(object)), "(0,1)"))
    }
    bf <- getBernsteinFunction(object)
    if (!test_class(bf@original, "AlphaStableBernsteinFunction")) {
      return(sprintf(ERR_MSG_BF_TYPE, classLabel(class(object))))
    }

    invisible(TRUE)
  })


#' @importFrom checkmate qassert
#' @importFrom rmo SumOfBernsteinFunctions AlphaStableBernsteinFunction
setMethod("constructBernsteinFunction", "AlphaStableExtMO2FParam",
  function(object, lambda, nu) {
    qassert(lambda, "N1(0,)")
    qassert(nu, "N1(0,1)")
    ScaledBernsteinFunction(
      scale = lambda,
      original = AlphaStableBernsteinFunction(alpha = nu)
    )
  })

#' @importFrom checkmate qassert
setMethod("calcAlpha2Nu", "AlphaStableExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,1)")
    log2(2 - value)
  })

#' @importFrom checkmate qassert
setMethod("calcNu2Alpha", "AlphaStableExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,1)")
    2 - 2 ^ value
  })

setReplaceMethod("setNu", "AlphaStableExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,1)")
    callNextMethod()
  })



#' @rdname ExtMO2FParam-class
#'
#' @section Poisson calibration parameter class:
#' Corresponds to a Lévy subrodinator which is a convex combination of a
#' Poisson subordinator with jump size `nu` and a pure-drift subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = \operatorname{e}^{-\nu}x + (1 - \operatorname{e}^{-x \nu})}
#'  \item \eqn{\nu = -log(1 - sqrt(\alpha))} and \eqn{\alpha = (1 - \operatorname{e}^{-\nu})^2}
#' }
#'
#' @export PoissonExtMO2FParam
PoissonExtMO2FParam <- setClass("PoissonExtMO2FParam", # nolint
  contains = "ExtMO2FParam")


#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions LinearBernsteinFunction
#'   PoissonBernsteinFunction
#' @importFrom checkmate qtest test_class
setValidity("PoissonExtMO2FParam",
  function(object) {
    if (!qtest(object@nu, "N1[0,)")) {
      return(sprintf(ERR_MSG_NU1_INTERVAL, classLabel(class(object)), "[0,)"))
    }
    bf <- getBernsteinFunction(object)
    if (!(
        test_class(bf@original, "SumOfBernsteinFunctions") &&
        test_class(bf@original@first, "LinearBernsteinFunction") &&
        test_class(bf@original@second, "PoissonBernsteinFunction"))) {
      return(sprintf(ERR_MSG_BF_TYPE, classLabel(class(object))))
    }

    invisible(TRUE)
  })


#' @importFrom rmo SumOfBernsteinFunctions LinearBernsteinFunction PoissonBernsteinFunction
#' @importFrom checkmate qassert
setMethod("constructBernsteinFunction", "PoissonExtMO2FParam",
  function(object, lambda, nu) {
    qassert(lambda, "N1(0,)")
    qassert(nu, "N1[0,)")
    ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = exp(-nu)),
        second = PoissonBernsteinFunction(lambda = 1, eta = nu)
      )
    )
  })

#' @importFrom checkmate qassert
setMethod("calcAlpha2Nu", "PoissonExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1)")
    -log(1 - sqrt(value))
  })

#' @importFrom checkmate qassert
setMethod("calcNu2Alpha", "PoissonExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,)")
    (1 - exp(-value))^2
  })

setReplaceMethod("setNu", "PoissonExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,)")
    callNextMethod()
  })



#' @rdname ExtMO2FParam-class
#'
#' @section Exponential calibration parameter class:
#' Corresponds to a Lévy subordinator which is a convex combination of an
#' Exponential-jump compound Poisson process with rate `nu` and unit-intensity
#' and a pure-drift subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = (1 - 1 / (1 + \nu))x + x   / (x + \nu)}
#'   \item \eqn{\nu = 0.5 \cdot (-3 + \sqrt{1 + 8 / \alpha})}
#'     and \eqn{\alpha = 2 / (1 + \nu) - 2 / (2 + \nu)}
#' }
#'
#' @export ExponentialExtMO2FParam
ExponentialExtMO2FParam <- setClass("ExponentialExtMO2FParam", # nolint
  contains = "ExtMO2FParam")

#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions LinearBernsteinFunction
#'   ExponentialBernsteinFunction
#' @importFrom checkmate qtest test_class
setValidity("ExponentialExtMO2FParam",
  function(object) {
    if (!qtest(object@nu, "N1(0,)")) {
      return(sprintf(ERR_MSG_NU1_INTERVAL, classLabel(class(object)), "(0,)"))
    }
    bf <- getBernsteinFunction(object)
    if (!(
        test_class(bf@original, "SumOfBernsteinFunctions") &&
        test_class(bf@original@first, "LinearBernsteinFunction") &&
        test_class(bf@original@second, "ExponentialBernsteinFunction"))) {
      return(sprintf(ERR_MSG_BF_TYPE, classLabel(class(object))))
    }

    invisible(TRUE)
  })


#' @importFrom rmo SumOfBernsteinFunctions LinearBernsteinFunction ExponentialBernsteinFunction
#' @importFrom checkmate qassert
setMethod("constructBernsteinFunction", "ExponentialExtMO2FParam",
  function(object, lambda, nu, ...) {
    qassert(lambda, "N1(0,)")
    qassert(nu, "N1(0,)")
    ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = 1 - 1 / (1 + nu)),
        second = ExponentialBernsteinFunction(lambda = nu)
      )
    )
  })

#' @importFrom checkmate qassert
setMethod("calcAlpha2Nu", "ExponentialExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,1)")
    0.5 * (-3 + sqrt(1 + 8 / value))
  })

#' @importFrom checkmate qassert
setMethod("calcNu2Alpha", "ExponentialExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    2 * (1 / (1 + value) - 1 / (2 + value))
  })

setReplaceMethod("setNu", "ExponentialExtMO2FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    callNextMethod()
  })
