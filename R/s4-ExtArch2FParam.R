#' @include s4-ExtMO2FParam.R checkmate.R
NULL

# nolint start
ERR_MSG_COPULA_TYPE <- "`copula` is of wrong type"
# nolint end

#' Two-factor extendible Archimedean calibration parameter classes
#'
#' @description
#' [CalibrationParam-class]-class with two parameters for the extendible
#' Archimedean-(survival-)copula model with Exponential margins for the
#' *(average) default counting process*.
#'
#' @slot lambda A non-negative number for the marginal rate.
#' @slot nu A numeric number for the model specific dependence parameter.
#'
#' @details
#' The model is defined by the assumption that the multivariate default times
#' \eqn{\tau = (\tau_1, \ldots, \tau_d)} are from an extendible
#' Archimedean-(survival-)copula with Exponential margins.
#' The (internal) dependency parameter \eqn{\nu} (model specific) has a
#' one-to-one relationship and can be replaced by *Spearman's Rho* `rho` (except
#' for the family based on Joe's copula) or *Kendall's Tau* `tau`. The possible
#' range for `rho` and `tau` is from zero to one (boundaries might not be
#' included).
#'
#' @importFrom copula archmCopula
#'
#' @export ExtArch2FParam
ExtArch2FParam <- setClass("ExtArch2FParam", # nolint
  contains = "CalibrationParam",
  slots = c(lambda = "numeric", nu = "numeric",
  family = "character", copula = "archmCopula", survival = "logical"))


#' @importFrom copula archmCopula
setMethod("getDimension", "archmCopula",
  function(object) {
    object@dimension
  })

#' @importFrom copula getTheta
setMethod("getNu", "archmCopula",
  function(object) {
    getTheta(object)
  })

setGeneric("getFamily",
  function(object) {
    standardGeneric("getFamily")
  })
setMethod("getFamily", "ExtArch2FParam",
  function(object) {
    object@family
  })

setGeneric("getCopula",
  function(object) {
    standardGeneric("getCopula")
  })
setMethod("getCopula", "ExtArch2FParam",
  function(object) {
    object@copula
  })

setGeneric("setCopula<-",
  function(object, value) {
    standardGeneric("setCopula<-")
  })
#' @importFrom checkmate assert_class
setReplaceMethod("setCopula", "ExtArch2FParam",
  function(object, value) {
    assert_class(value, paste0(tolower(getFamily(object)), "Copula"))
    object@copula <- value
    if (!isTRUE(getDimension(object) == getDimension(value)) && !is.na(getNu(value))) {
      setDimension(object) <- getDimension(value)
    }
    if (!isTRUE(getNu(object) == getNu(value)) && !is.na(getNu(value))) {
      setNu(object) <- getNu(value)
    }

    invisible(object)
  })

setGeneric("setFamily<-",
  function(object, value) {
    standardGeneric("setFamily<-")
  })
setReplaceMethod("setFamily", "ExtArch2FParam",
  function(object, value) {
    assert_choice(value, c("Clayton", "Frank", "Gumbel", "Joe"))
    object@family <- value
    setCopula(object) <- constructCopula(object)

    invisible(object)
  })

setGeneric("constructCopula",
  function(object, ...) {
    standardGeneric("constructCopula")
  })
#' @importFrom copula archmCopula
#' @importFrom checkmate assert_choice test_class qtest qassert
setMethod("constructCopula", "ExtArch2FParam",
  function(object, family = getFamily(object), nu = getNu(object), d = getDimension(object)) {
    assert_choice(family, c("Clayton", "Frank", "Gumbel", "Joe"))
    out <- getCopula(object)
    if (test_class(out, paste0(tolower(family), "Copula")) && isTRUE(getDimension(out) == d)) {
      out <- setTheta(out, nu)
    } else if (missing(nu) && missing(d)) {
      if (qtest(nu, "N1") && qtest(d, "X1[2,)")) {
        out <- archmCopula(family, nu, d)
      } else {
        out <- archmCopula(family)
      }
    } else {
      qassert(nu, "N1")
      qassert(d, "X1[2,)")
      out <- archmCopula(family, nu, d)
    }

    out
  })

#' @importFrom copula archmCopula
#' @importFrom checkmate qassert qtest test_choice
setReplaceMethod("setDimension", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "X1[2,)")
    object <- callNextMethod()
    setCopula(object) <- constructCopula(object)

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
#' @importFrom copula setTheta archmCopula
#' @importFrom checkmate qassert qtest test_choice
setReplaceMethod("setNu", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1")
    object@nu <- value
    setCopula(object) <- constructCopula(object)

    invisible(object)
  })

#' @importFrom copula rho
setMethod("getRho", "ExtArch2FParam",
  function(object) {
    rho(getCopula(object))
  })
#' @importFrom checkmate qassert
setReplaceMethod("setRho", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invRho(object, value)

    invisible(object)
  })

#' @importFrom copula tau
setMethod("getTau", "ExtArch2FParam",
  function(object) {
    tau(getCopula(object))
  })
#' @importFrom checkmate qassert
setReplaceMethod("setTau", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invTau(object, value)

    invisible(object)
  })

setGeneric("getSurvival",
  function(object) {
    standardGeneric("getSurvival")
  })
setMethod("getSurvival", "ExtArch2FParam",
  function(object) {
    object@survival
  })

setGeneric("setSurvival<-",
  function(object, value) {
    standardGeneric("setSurvival<-")
  })
#' @importFrom checkmate qassert
setReplaceMethod("setSurvival", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "B1")
    object@survival <- value

    invisible(object)
  })


#' @importFrom copula getTheta
#' @importFrom checkmate qtest test_class test_choice
setValidity("ExtArch2FParam", # nolint
  function(object) {
    if (!qtest(object@lambda, "N1(0,)")) {
      return(sprintf(ERR_MSG_LAMBDA, classLabel(class(object))))
    }
    if (!qtest(object@nu, "N1")) {
      return(sprintf(ERR_MSG_NU1_INTERVAL, classLabel(class(object)), "[0,1]"))
    }
    if (!(
        qtest(object@survival, "B1") &&
        test_class(object@copula, "archmCopula") &&
        test_choice(object@family, c("Clayton", "Frank", "Gumbel", "Joe")) &&
        test_choice(getDimension(object@copula), getDimension(object)) &&
        test_choice(getNu(object@copula), getNu(object)))) {
      return(ERR_MSG_COPULA_TYPE)
    }

    invisible(TRUE)
  })


#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,ExtArch2FParam-method
#' @aliases initialize,ExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param dim Dimension.
#' @param lambda Marginal intensity.
#' @param nu (Internal) dependence parameter (see [copula::archmCopula-class]).
#' @param rho Bivariate Spearman's Rho.
#' @param tau Bivariate Kendall's Tau.
#' @param survival Flag if survival copula is specified (default, except for Clayton).
#' @param family Name of the Archimedean copula family (see [copula::archmCopula-class]).
#'
#' @examples
#' ExtArch2FParam(dim = 5L, lambda = 8e-2, rho = 4e-1, family = "Clayton", survival = TRUE)
#' ExtArch2FParam(dim = 5L, lambda = 8e-2, rho = 4e-1, family = "Gumbel", survival = FALSE)
#' ExtArch2FParam(dim = 5L, lambda = 8e-2, tau = 4e-1, family = "Frank", survival = FALSE)
#' ExtArch2FParam(dim = 5L, lambda = 8e-2, tau = 4e-1, family = "Joe", survival = FALSE)
#' @importFrom copula archmCopula
setMethod("initialize", "ExtArch2FParam",
  function(.Object, # nolint
      dim, lambda, nu, rho, tau, survival,
      family = c("Clayton", "Frank", "Gumbel", "Joe")) {
    setFamily(.Object) <- match.arg(family)
    setSurvival(.Object) <- survival
    if (!missing(dim) && !missing(lambda) &&
          !(missing(nu) && missing(rho) && missing(tau))) {
      setDimension(.Object) <- dim
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


#' @importFrom copula iRho
#' @importFrom checkmate qassert
setMethod("invRho", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    iRho(getCopula(object), value)
  })

#' @importFrom copula iTau
#' @importFrom checkmate qassert
setMethod("invTau", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    iTau(getCopula(object), value)
  })



#' @describeIn ExtArch2FParam-class
#'    simulates the vector of *default times* and returns a matrix `x` with
#'    `dim(x) == c(n_sim, getDimension(object))`.
#' @aliases simulate_dt,ExtArch2FParam-method
#'
#' @inheritParams simulate_dt
#' @param n_sim Number of samples.
#'
#' @section Simulation:
#' The default times are sampled in a two-stage procedure: First a sample is
#' drawn from the Archimedean copula, see [copula::archmCopula-class] and
#' [copula::rCopula()]; then the results are transformed using [stats::qexp()].
#'
#' @examples
#' parm <- FrankExtArch2FParam(dim = 5L, lambda = 8e-2, rho = 4e-1)
#' simulate_dt(parm, n_sim = 5L)
#'
#' @importFrom stats qexp
#' @importFrom copula rCopula
#' @include utils.R
setMethod("simulate_dt", "ExtArch2FParam",
  function(object, ..., n_sim = 1e4L) {
    out <- qexp(
      rCopula(n_sim, getCopula(object)),
      rate = getLambda(object), lower.tail = !getSurvival(object))

    simplify2vector(out)
  })

#' @describeIn ExtArch2FParam-class
#'   calculates the *payoff equation* for a *portfolio CDS* (vectorized w.r.t.
#'   the argumentes `recovery_rate`, `coupon`, and `upfront`).
#' @aliases expected_pcds_equation,ExtArch2FParam-method
#'
#' @inheritParams expected_pcds_equation
#'
#' @inheritSection ExtMO2FParam-class Expected portfolio CDS loss
#'
#' @examples
#' parm <- FrankExtArch2FParam(75L, 8e-2, rho = 4e-1)
#' expected_pcds_equation(
#'   parm, times =seq(0, 5, by = 0.25), discount_factors = rep(1, 21L), recovery_rate = 0.4,
#'   coupon = 1e-1, upfront = 0)
#' expected_pcds_equation(
#'   parm, times =seq(0, 5, by = 0.25), discount_factors = rep(1, 21L), recovery_rate = 0.4,
#'   coupon = 1e-1, upfront = 0, method = "mc", n_sim = 1e1)
#'
#' @importFrom stats pexp
#' @importFrom checkmate qassert assert check_numeric
#' @importFrom vctrs vec_size_common vec_recycle
#' @include RcppExports.R
#'
#' @export
setMethod("expected_pcds_equation", "ExtArch2FParam",
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


#' @describeIn ExtArch2FParam-class Display the object.
#' @aliases show,ExtArch2FParam-method
#'
#' @inheritParams methods::show
#'
#' @export
setMethod("show", "ExtArch2FParam",
 function(object) {
   cat(sprintf("An object of class %s\n", classLabel(class(object))))
   if (isTRUE(validObject(object, test = TRUE))) {
     cat(sprintf("Dimension: %i\n", getDimension(object)))
     cat("Parameter:\n")
     cat(sprintf("* %s: %s\n", "Lambda", format(getLambda(object))))
     if (isTRUE(getFamily(object) != "Joe")) {
       cat(sprintf("* %s: %s\n", "Rho", format(getRho(object))))
     }
     cat(sprintf("* %s: %s\n", "Tau", format(getTau(object))))
     cat("Internal parameter:\n")
     cat(sprintf("* %s: %s\n", "Nu", format(getNu(object))))
     cat(sprintf("* %s: %s\n", "Survival copula", ifelse(getSurvival(object), "Yes", "No")))
   } else {
     cat("\t (invalid or not initialized)\n")
   }

   invisible(NULL)
 })




#' @rdname ExtArch2FParam-class
#'
#' @export ClaytonExtArch2FParam
ClaytonExtArch2FParam <- setClass("ClaytonExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "claytonCopula"))

#' @importFrom checkmate assert_choice
setValidity("ClaytonExtArch2FParam",
  function(object) {
    assert_choice(object@family, "Clayton")

    invisible(TRUE)
  })

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,ClaytonExtArch2FParam-method
#' @aliases initialize,ClaytonExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' ClaytonExtArch2FParam()
#' ClaytonExtArch2FParam(5L, 8e-2, rho = 4e-1)
#' ClaytonExtArch2FParam(5L, 8e-2, tau = 4e-1)
setMethod("initialize", "ClaytonExtArch2FParam",
  function(.Object, ..., survival = FALSE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Clayton"))
  })

#' @importFrom checkmate assert_choice
setReplaceMethod("setFamily", "ClaytonExtArch2FParam",
  function(object, value) {
    assert_choice(value, "Clayton")
    object <- callNextMethod()

    invisible(object)
  })



#' @rdname ExtArch2FParam-class
#'
#' @export FrankExtArch2FParam
FrankExtArch2FParam <- setClass("FrankExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "frankCopula"))

#' @importFrom checkmate assert_choice
setValidity("FrankExtArch2FParam",
  function(object) {
    assert_choice(object@family, "Frank")

    invisible(TRUE)
  })

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,FrankExtArch2FParam-method
#' @aliases initialize,FrankExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' FrankExtArch2FParam()
#' FrankExtArch2FParam(5L, 8e-2, rho = 4e-1)
#' FrankExtArch2FParam(5L, 8e-2, tau = 4e-1)
setMethod("initialize", "FrankExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Frank"))
  })

#' @importFrom checkmate assert_choice
setReplaceMethod("setFamily", "FrankExtArch2FParam",
  function(object, value) {
    assert_choice(value, "Frank")
    object <- callNextMethod()

    invisible(object)
  })



#' @rdname ExtArch2FParam-class
#'
#' @export GumbelExtArch2FParam
GumbelExtArch2FParam <- setClass("GumbelExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "gumbelCopula"))

#' @importFrom checkmate assert_choice
setValidity("GumbelExtArch2FParam",
  function(object) {
    assert_choice(object@family, "Gumbel")

    invisible(TRUE)
  })

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,GumbelExtArch2FParam-method
#' @aliases initialize,GumbelExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' GumbelExtArch2FParam()
#' GumbelExtArch2FParam(5L, 8e-2, rho = 4e-1)
#' GumbelExtArch2FParam(5L, 8e-2, tau = 4e-1)
setMethod("initialize", "GumbelExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Gumbel"))
  })

#' @importFrom checkmate assert_choice
setReplaceMethod("setFamily", "GumbelExtArch2FParam",
  function(object, value) {
    assert_choice(value, "Gumbel")
    object <- callNextMethod()

    invisible(object)
  })



#' @rdname ExtArch2FParam-class
#'
#' @export JoeExtArch2FParam
JoeExtArch2FParam <- setClass("JoeExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "joeCopula"))

#' @importFrom checkmate assert_choice
setValidity("JoeExtArch2FParam",
  function(object) {
    assert_choice(object@family, "Joe")

    invisible(TRUE)
  })

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,JoeExtArch2FParam-method
#' @aliases initialize,JoeExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' JoeExtArch2FParam()
#' JoeExtArch2FParam(5L, 8e-2, tau = 4e-1)
setMethod("initialize", "JoeExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Joe"))
  })

#' @importFrom checkmate assert_choice
setReplaceMethod("setFamily", "JoeExtArch2FParam",
  function(object, value) {
    assert_choice(value, "Joe")
    object <- callNextMethod()

    invisible(object)
  })

setReplaceMethod("setRho", "JoeExtArch2FParam",
  function(object, value) {
    stop(paste0("Spearman's Rho not implemented for family ",
      getFamily(object), "ExtArch2FParam"))
  })
