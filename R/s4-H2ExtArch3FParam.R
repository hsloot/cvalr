#' @include s4-H2ExCalibrationParam.R s4-H2ExtMO3FParam.R s4-H2ExtGaussian3FParam.R checkmate.R
NULL

#' Three-factor H2-extendible Archimedean calibration parameter
#'
#' [CalibrationParam-class] for the H2-extendible Archimedean copula with Exponential margin model
#' for the *(average) default counting process*  with 3 parameter. Extends
#' [H2ExCalibrationParam-class] and related to [ExtArch2FParam-class].
#'
#' @slot lambda A non-negative number for the marginal rate.
#' @slot nu A numeric vector of length 2 for the model specific dependence parameters (global and
#'   component specific; range depends on specific model). Use `rho` or `tau` in the constructor to
#'   set dependence parameter.
#'
#' @details
#' The model is defined by the assumption that the *multivariate default times* \eqn{\tau = (\tau_1,
#' \ldots, \tau_d)} are from a H2-extendible Archimedean copula model with Exponential margins. Per
#' default, the Archimedean copula is used as a survival copula, except for the Clayton-family.
#' The model is specified by three parameters (in addition to the composition): The *marginal rate*
#' `lambda` and the (internal) *outer* and *inner dependency parameters* `nu` (see
#' [outer_nacopula-class]). The dependency parameter `nu` should not be set by the uesr; instead
#' they should provide either `rho` (*Spearman's Rho*) or `tau` (*Kendall's Tau*).
#' The parameters `rho` or `tau` should be between zero and one, of length 2, and non-decreasing;
#' the first value represents the *outer dependence* between components of different partition
#' elements and the second value represents the *inner dependence* between components of the same
#' partition element.
#'
#' For details on the underlying extendible model, see [ExtArch2FParam-class].
#'
#' @export H2ExtArch3FParam
H2ExtArch3FParam <- setClass("H2ExtArch3FParam", # nolint
  contains = c("H2ExCalibrationParam"),
  slots = c(lambda = "numeric", nu = "numeric", family = "character",
    survival = "logical", copula = "outer_nacopula"))


setMethod("getModelName", "H2ExtArch3FParam",
  function(object) {
    "ExtArch2FParam"
  })

setMethod("getFamily", "outer_nacopula",
  function(object) {
    object@copula@name
  })

setMethod("getDimension", "outer_nacopula",
  function(object) {
    dim(object)
  })

setMethod("getComposition", "outer_nacopula",
  function(object) {
    map_int(object@childCops, dim)
  })

#' @importFrom copula getTheta
setMethod("getNu", "outer_nacopula",
  function(object) {
    c(getTheta(object@copula), getTheta(object@childCops[[1]]@copula))
  })

setMethod("getFamily", "H2ExtArch3FParam",
  function(object) {
    object@family
  })

setMethod("getSurvival", "H2ExtArch3FParam",
  function(object) {
    object@survival
  })

#' @importFrom checkmate qassert
setReplaceMethod("setSurvival", "H2ExtArch3FParam",
  function(object, value) {
    qassert(value, "B1")
    object@survival <- value

    invisible(object)
  })

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

setMethod("getCopula", "H2ExtArch3FParam",
  function(object) {
    object@copula
  })

#' @importFrom checkmate assert check_class check_choice
setReplaceMethod("setCopula", "H2ExtArch3FParam", {
  function(object, value) {
    assert(combine = "and",
      check_class(value, "outer_nacopula"), check_choice(value@copula@name, getFamily(object)))
    object@copula <- value

    invisible(object)
  }
})

#' @importFrom copula onacopulaL
#' @importFrom purrr map
#' @importFrom checkmate test_numeric qtest assert_numeric qassert
setMethod("constructCopula", "H2ExtArch3FParam",
  function(object,
           family = getFamily(object), nu = getNu(object), composition = getComposition(object)) {
    assert_choice(family, c("Clayton", "Frank", "Gumbel", "Joe"))
    out <- getCopula(object)
    if (!missing(family) || !missing(nu) || !missing(composition)) {
      assert_numeric(nu, any.missing = FALSE, finite = TRUE, len = 2L, sorted = TRUE)
      qassert(composition, "X+[1,)")
      qassert(sum(composition), "X1[2,)")
      nac_list <- list(nu[[1]], c(), map(getPartition(object), ~list(nu[[2]], .)))
      out <- onacopulaL(family, nac_list)
    } else if (
        test_numeric(nu, any.missing = FALSE, finite = TRUE, len = 2L, sorted = TRUE) &&
        qtest(composition, "X+[1,)") && qtest(sum(composition), "X1[2,)")) {
      nac_list <- list(nu[[1]], c(), map(getPartition(object), ~list(nu[[2]], .)))
      out <- onacopulaL(family, nac_list)
    }

    out
  })

setReplaceMethod("setFamily", "H2ExtArch3FParam",
  function(object, value) {
    assert_choice(value, c("Clayton", "Frank", "Gumbel", "Joe"))
    object@family <- value
    copula <- constructCopula(object)
    if (
        test_class(copula, "outer_nacopula") &&
        test_choice(copula@copula@name, getFamily(object))) {
      setCopula(object) <- copula
    }
    invisible(object)
  })

setReplaceMethod("setComposition", "H2ExtArch3FParam",
  function(object, value) {
    qassert(value, "X+[1,)")
    dim <- sum(value)
    qassert(dim, "X1[2,)")
    object <- callNextMethod()
    copula <- constructCopula(object)
    if (
        test_class(copula, "outer_nacopula") &&
        test_choice(copula@copula@name, getFamily(object))) {
      setCopula(object) <- copula
    }

    invisible(object)
  })

setMethod("getNu", "H2ExtArch3FParam",
  function(object) {
    object@nu
  })

#' @importFrom purrr map
#' @importFrom checkmate qassert
setReplaceMethod("setNu", "H2ExtArch3FParam",
  function(object, value) {
    qassert(value, "N2")
    object@nu <- value
    copula <- constructCopula(object)
    if (
        test_class(copula, "outer_nacopula") &&
        test_choice(copula@copula@name, getFamily(object))) {
      setCopula(object) <- copula
    }

    invisible(object)
  })

#' @importFrom copula rho
#' @importFrom purrr map
#' @importFrom checkmate assert_numeric
setMethod("calcNu2Rho", "H2ExtArch3FParam",
  function(object, value) {
    assert_numeric(value, any.missing = FALSE, len = 2L, sorted = TRUE)
    map_dbl(value, ~rho(archmCopula(getFamily(object), .)))
  })

#' @importFrom copula rho
#' @importFrom purrr map_dbl
setMethod("getRho", "H2ExtArch3FParam",
  function(object) {
    calcNu2Rho(object, getNu(object))
  })

#' @importFrom copula iRho
#' @importFrom purrr map
#' @importFrom checkmate assert_numeric
setMethod("calcRho2Nu", "H2ExtArch3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    map_dbl(value, ~iRho(archmCopula(getFamily(object)), .))
  })

#' @importFrom checkmate assert_numeric
setMethod("invRho", "H2ExtArch3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    calcRho2Nu(object, value)
  })

#' @importFrom checkmate assert_numeric
setReplaceMethod("setRho", "H2ExtArch3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    setNu(object) <- invRho(object, value)

    invisible(object)
  })

#' @importFrom copula tau
#' @importFrom purrr map
#' @importFrom checkmate assert_numeric
setMethod("calcNu2Tau", "H2ExtArch3FParam",
  function(object, value) {
    assert_numeric(value, any.missing = FALSE, len = 2L, sorted = TRUE)
    map_dbl(value, ~tau(archmCopula(getFamily(object), .)))
  })

#' @importFrom copula rho
#' @importFrom purrr map_dbl
setMethod("getTau", "H2ExtArch3FParam",
  function(object) {
    calcNu2Tau(object, getNu(object))
  })

#' @importFrom copula iTau
#' @importFrom purrr map
#' @importFrom checkmate assert_numeric
setMethod("calcTau2Nu", "H2ExtArch3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    map_dbl(value, ~iTau(archmCopula(getFamily(object)), .))
  })

#' @importFrom checkmate assert_numeric
setMethod("invTau", "H2ExtArch3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    calcTau2Nu(object, value)
  })

#' @importFrom checkmate assert_numeric
setReplaceMethod("setTau", "H2ExtArch3FParam",
  function(object, value) {
    assert_numeric(value, lower = 0, upper = 1, any.missing = FALSE, len = 2L, sorted = TRUE)
    setNu(object) <- invTau(object, value)

    invisible(object)
  })


#' @importFrom checkmate qtest test_class test_choice
setValidity("H2ExtArch3FParam", # nolint
  function(object) {
    if (!qtest(object@lambda, "N1(0,)")) {
      return(ERR_MSG_LAMBDA)
    }
    if (!qtest(object@nu, "N2(0,)")) {
      return(sprintf(ERR_MSG_NU2_INTERVAL, "(0,)"))
    }
    if (!(
        qtest(object@survival, "B1") &&
        test_class(object@copula, "outer_nacopula") &&
        test_choice(object@family, c("Clayton", "Frank", "Gumbel", "Joe")) &&
        test_choice(getDimension(object@copula), getDimension(object)) &&
        test_equal(getNu(object@copula), getNu(object)))) {
      return(ERR_MSG_COPULA_TYPE)
    }

    invisible(TRUE)
  })


#' @describeIn H2ExtArch3FParam-class Constructor
#' @aliases initialize,H2ExtArch3FParam-method
#' @aliases initialize,H2ExtArch3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param composition Positive integerish vector for the component-composition.
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
#'   composition = c(3L, 3L, 4L, 5L), lambda = 8e-2, tau = c(3e-1, 5e-1), survival = TRUE)
#'
#' @importFrom utils head
#' @importFrom copula archmCopula onacopulaL iRho iTau
#' @importFrom purrr imap map map2
setMethod("initialize", "H2ExtArch3FParam",
  function(.Object, # nolint
      composition = c(2L, 3L),
      lambda = 1e-1, nu = c(0.2, 0.3), rho = NULL, tau = NULL,
      survival = TRUE,
      family = c("Clayton", "Frank", "Gumbel", "Joe")) {
    setFamily(.Object) <- match.arg(family)
    setSurvival(.Object) <- survival
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


#' @describeIn H2ExtArch3FParam-class
#'   returns the expected portfolio CDS loss for a specific time-point.
#' @aliases expected_pcds_loss,H2ExtArch3FParam-method
#'
#' @inheritParams probability_distribution
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#'
#' @inheritSection ExtMO2FParam-class Expected portfolio CDS loss
#'
#' @examples
#' parm <- ClaytonH2ExtArch3FParam(c(3, 3, 4, 5), 8e-2, rho = c(3e-1, 6e-1))
#' expected_pcds_loss(parm, times = 0.25, recovery_rate = 0.4)
#' expected_pcds_loss(parm, times = seq(0, 1, by = 0.25), recovery_rate = 0.4)
#' expected_pcds_loss(parm, times = seq(0, 1, by = 0.25), recovery_rate = 0.4,
#'   method = "CalibrationParam")
#' expected_pcds_loss(parm, times = seq(0, 1, by = 0.25), recovery_rate = 0.4,
#'   method = "CalibrationParam",
#'   pd_args = list(method = "CalibrationParam", seed = 1623,
#'     sim_args = list(n_sim = 1e2L)))
#'
#' @importFrom stats pexp
#' @importFrom checkmate qassert
#'
#' @export
setMethod("expected_pcds_loss", "H2ExtArch3FParam",
  function(object, times, recovery_rate, ...,
      method = c("default", "H2ExtArch3FParam", "CalibrationParam")) {
    method <- match.arg(method)
    if (isTRUE("default" == method || "H2ExtArch3FParam" == method)) {
      qassert(times, "N+[0,)")
      qassert(recovery_rate, "N1[0,1]")
      out <- (1 - recovery_rate) * pexp(times, rate = getLambda(object))
    } else {
      out <- callNextMethod(object, times, recovery_rate, ..., method = method)
    }

    out
  })


#' @describeIn H2ExtArch3FParam-class Display the object.
#' @aliases show,H2ExtArch3FParam-method
#'
#' @inheritParams methods::show
#'
#' @export
setMethod("show", "H2ExtArch3FParam",
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
     if (isTRUE("Joe" != getFamily(object))) {
       cat(sprintf("* %s: %s\n", "Rho", to_vector(format(getRho(object)))))
     }
     cat(sprintf("* %s: %s\n", "Tau", to_vector(format(getTau(object)))))
     cat("Internal parameter:\n")
     cat(sprintf("* %s: %s\n", "Nu", to_vector(format(getNu(object)))))
   } else {
     cat("\t (invalid or not initialized)\n")
   }

   invisible(NULL)
  })


#' @describeIn H2ExtArch3FParam-class
#'    simulates the vector of *default times* and returns a matrix `x` with
#'    `dim(x) == c(n_sim, getDimension(object))`.
#' @aliases simulate_dt,H2ExtArch3FParam-method
#'
#' @inheritParams simulate_dt
#' @param n_sim Number of samples.
#'
#' @section Simulation:
#' The default times are sampled in a two-stage procedure: First a sample is drawn from the
#' [copula::outer_nacopula-class] copula whose dependence reflect the inner- and outer-dependency
#' parameters; then the results are transformed using [stats::qexp()].
#'
#' @examples
#' parm <- FrankH2ExtArch3FParam(composition = c(2L, 4L, 2L), lambda = 8e-2, rho = c(2e-1, 7e-1))
#' simulate_dt(parm, n_sim = 5L)
#'
#' @importFrom stats qexp
#' @importFrom copula rCopula
#' @include utils.R
setMethod("simulate_dt", "H2ExtArch3FParam",
  function(object, ..., n_sim = 1e4L) {
    rCopula(n_sim, getCopula(object)) %>%
      qexp(rate = getLambda(object), lower.tail = !getSurvival(object))
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
#' ClaytonH2ExtArch3FParam()
#' ClaytonH2ExtArch3FParam(composition = c(3L, 3L, 4L, 5L), lambda = 8e-2, tau = c(3e-1, 5e-1))
#' ClaytonH2ExtArch3FParam(composition = c(3L, 3L, 4L, 5L), lambda = 8e-2, rho = c(3e-1, 5e-1))
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
#' FrankH2ExtArch3FParam()
#' FrankH2ExtArch3FParam(composition = c(3L, 3L, 4L, 5L), lambda = 8e-2, tau = c(3e-1, 5e-1))
#' FrankH2ExtArch3FParam(composition = c(3L, 3L, 4L, 5L), lambda = 8e-2, rho = c(3e-1, 5e-1))
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
#' GumbelH2ExtArch3FParam()
#' GumbelH2ExtArch3FParam(composition = c(3L, 3L, 4L, 5L), lambda = 8e-2, tau = c(3e-1, 5e-1))
#' GumbelH2ExtArch3FParam(composition = c(3L, 3L, 4L, 5L), lambda = 8e-2, rho = c(3e-1, 5e-1))
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
#' JoeH2ExtArch3FParam()
#' JoeH2ExtArch3FParam(composition = c(3L, 3L, 4L, 5L), lambda = 8e-2, tau = c(3e-1, 5e-1))
setMethod("initialize", "JoeH2ExtArch3FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Joe"))
  })


setMethod("getModelName", "JoeH2ExtArch3FParam",
  function(object) {
    "JoeExtArch2FParam"
  })

setReplaceMethod("setRho", "JoeH2ExtArch3FParam",
  function(object, value) {
    stop(paste0("Spearman's Rho not implemented for family ",
      getFamily(object), "H2ExtArch3FParam"))
  })
