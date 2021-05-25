#' @include s4-ExMarkovParam.R checkmate.R
NULL

# nolint start
ERR_MSG_EXINTENSITIES <- "`ex_intensities` must be (scaled) exchangeable MO intensity vector"
# nolint end

#' Exchangeable Marshall--Olkin calibration parameter
#'
#' @description
#' [CalibrationParam-class]-class for the exchangeable Marshall-Olkin model for
#' the *average default counting process*. Extends [ExMarkovParam-class].
#'
#' @slot ex_intensities The (scaled) exchangeable intensities  of the
#'   exchangeable Marshall-Olkin distribution (see [rmo::exIntensities()]).
#'
#' @details
#' The model is defined by the assumption that the multivariate default times
#' \eqn{\tau = (\tau_1, \ldots, \tau_d)} are exchangeable Marshall-Olkin.
#' The joint survival function of all portfolio items is
#' \deqn{
#'   P(\tau > t)
#'     = \exp{(- a_0 t_{[1]} - \cdots - a_{d-1} t_{[d]})} ,
#' }
#' for \eqn{t_{[1]} \geq \cdots \geq t_{[d]}} begin the descending version of
#' \eqn{t} and
#' \deqn{
#'   a_{i}
#'     = \sum_{l=0}^{d-i-1} \binom{d-i-1}{l} \lambda_{l+1} .
#' }
#' The (scaled) exchangeable intensities, provided to the constructor are
#' \deqn{
#'   \binom{d}{i} \lambda_{i} , \ i \in \{ 1 , \ldots , d \} .
#' }
#'
#' @export ExMOParam
ExMOParam <- setClass("ExMOParam", # nolint
  contains = "ExMarkovParam",
  slots = c(ex_intensities = "numeric"))


setGeneric("getExIntensities",
  function(object) {
    standardGeneric("getExIntensities")
  })
setMethod("getExIntensities", "ExMOParam",
  function(object) {
    object@ex_intensities
  })

setGeneric("setExIntensities<-",
  function(object, value) {
    standardGeneric("setExIntensities<-")
  })
#' @importFrom checkmate qassert assert_numeric
setReplaceMethod("setExIntensities", "ExMOParam",
  function(object, value) {
    assert_numeric(value, lower = 0, finite = TRUE, any.missing = FALSE, min.len = 2)
    qassert(max(value), "N1(0,)")
    object@ex_intensities <- value
    setExQMatrix(object) <- rmo:::exi2exqm(value)

    invisible(object)
  })


#' @importFrom checkmate qtest test_numeric
setValidity("ExMOParam",
  function(object) {
    if (!(test_numeric(
          object@ex_intensities, lower = 0, finite = TRUE, any.missing = FALSE,
          len = getDimension(object)) &&
        qtest(max(object@ex_intensities), "N1(0,)"))) {
      return(ERR_MSG_EXINTENSITIES)
    }

    invisible(TRUE)
  })


#' @describeIn ExMOParam-class Constructor
#' @aliases initialize,ExMOParam-method
#' @aliases initialize,ExMOParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ex_intensities The (scaled) exchangeable intensities  of the
#'   Marshall-Olkin distribution (see [rmo::exIntensities()]).
#'
#' @examples
#' ExMOParam()
#' ExMOParam(rmo::exIntensities(rmo::AlphaStableBernsteinFunction(0.4), 5L))
setMethod("initialize", "ExMOParam",
  definition = function(.Object, ex_intensities) { # nolint
    if (!missing(ex_intensities)) {
      setExIntensities(.Object) <- ex_intensities
      validObject(.Object)
    }

    invisible(.Object)
  })


#' @describeIn ExMOParam-class
#'    simulates the vector of *default times* and returns a matrix `x` with
#'    `dim(x) == c(n_sim, getDimension(object))`.
#' @aliases simulate_dt,ExMOParam-method
#'
#' @inheritParams simulate_dt
#' @param method Simulation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param n_sim Number of samples.
#'
#' @section Simulation:
#' The default times are sampled using [rmo::rexmo_markovian()].
#'
#'
#' @examples
#' parm <- ExMOParam(rmo::exIntensities(rmo::AlphaStableBernsteinFunction(0.4), 5L))
#' simulate_dt(parm, n_sim = 5L)
#'
#' @importFrom rmo rexmo_markovian
#' @include utils.R
#'
#' @export
setMethod("simulate_dt", "ExMOParam",
  function(object, ...,
      method = c("default", "ExMOParam", "ExMarkovParam"), n_sim = 1e4L) {
    method <- match.arg(method)
    if (isTRUE("default" == method)) {
      method <- "ExMOParam"
    }

    if (isTRUE("ExMOParam" == method)) {
      out <- rexmo_markovian(n_sim, getDimension(object), getExIntensities(object))
    } else {
      out <- callNextMethod(object, ..., method = method, n_sim = n_sim)
    }

    out
  })


#' @describeIn ExMOParam-class Display the object.
#' @aliases show,ExMOParam-method
#'
#' @export
setMethod("show", "ExMOParam",
  function(object) {
    cat(sprintf("An object of class %s\n", classLabel(class(object))))
    if (isTRUE(validObject(object, test = TRUE))) {
      cat(sprintf("Dimension: %i\n", getDimension(object)))
      cat("(Scaled) intensity vector:\n")
      capture.output(print(getExIntensities(object))) %>%
        paste0("\t", .) %>%
        writeLines
    } else {
      cat("\t (invalid or not initialized)\n")
    }


    invisible(NULL)
  })
