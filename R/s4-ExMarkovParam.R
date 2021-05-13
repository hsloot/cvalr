#' @include s4-CalibrationParam.R checkmate.R
NULL

#' Exchangeable Markovian calibration parameter
#'
#' @description
#' [CalibrationParam-class]-class for the exchangeable Markovian *default
#' counting process* model.
#'
#' @slot ex_qmatrix The \eqn{(d+1) \times (d+1)} Markov generator matrix of the
#'   default counting process.
#'
#' @export ExMarkovParam
ExMarkovParam <- setClass("ExMarkovParam", # nolint
  contains = "CalibrationParam",
  slots = c(ex_qmatrix = "matrix"))


setGeneric("getExQMatrix",
  function(object) {
    standardGeneric("getExQMatrix")
  })
setMethod("getExQMatrix", "ExMarkovParam",
  function(object) {
    object@ex_qmatrix
  })

setGeneric("setExQMatrix<-",
  function(object, value) {
    standardGeneric("setExQMatrix<-")
  })
#' @include checkmate.R
setReplaceMethod("setExQMatrix", "ExMarkovParam",
  function(object, value) {
    assert_exqmatrix(value, min.rows = 1L, min.cols = 1L)

    setDimension(object) <- nrow(value) - 1L
    object@ex_qmatrix <- value

    invisible(object)
  })


#' @include checkmate.R
setValidity("ExMarkovParam",
  function(object) {
    assert_exqmatrix(object@ex_qmatrix, nrows = object@dim+1L, ncols = object@dim + 1L)

    invisible(TRUE)
  })


#' @describeIn ExMarkovParam-class Constructor
#' @aliases initialize,ExMarkovParam-method
#' @aliases initialize,ExMarkovParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ex_qmatrix (Exchangeable) Generator-matrix of the
#'   *(average) default counting process*.
#'
#' @examples
#' ExMarkovParam(rmo::exQMatrix(rmo::AlphaStableBernsteinFunction(0.4), 5L))
setMethod("initialize", "ExMarkovParam",
  function(.Object, ex_qmatrix) { # nolint
    if (!missing(ex_qmatrix)) {
      setExQMatrix(.Object) <- ex_qmatrix
      validObject(.Object)
    }

    invisible(.Object)
  })


#' @describeIn ExMarkovParam-class
#'    simulates the vector of *default times* and returns a matrix `x` with
#'    `dim(x) == c(n_sim, getDimension(object))`.
#' @aliases simulate_dt,ExMarkovParam-method
#'
#' @inheritParams simulate_dt
#' @param method Simulation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param n_sim Number of samples.
#'
#' @examples
#' parm <- ExMarkovParam(rmo::exQMatrix(rmo::AlphaStableBernsteinFunction(0.4), 5L))
#' simulate_dt(parm, n_sim = 5L)
#'
#' @importFrom stats rexp
#' @include utils.R
#'
#' @export
setMethod("simulate_dt", "ExMarkovParam",
  function(object, ...,
      method = c("default", "ExMarkovParam"), n_sim = 1e1L) {
    method <- match.arg(method)
    out <- matrix(nrow = n_sim, ncol = object@dim)
    for (k in 1:n_sim) {
      state <- 0
      time <- 0
      while (state != object@dim) {
        wt <- rexp(1, rate = -object@ex_qmatrix[1+state, 1+state])
        time <- time + wt
        out[k, (1+state):object@dim] <- time
        state <- state +
          sample.int(
            n = object@dim-state, size = 1, replace = FALSE,
            prob = object@ex_qmatrix[1+state, (2+state):(object@dim+1)])
      }
      perm <- sample.int(n = object@dim, size = object@dim, replace = FALSE)
      out[k, perm] <- out[k, perm]
    }

    out
  })

#' @describeIn ExMarkovParam-class
#'   calculates the *probability vector* for the *average default count process*
#'   and returns a matrix `x` with `dim(x) == c(getDimension(object)+1L, length(times))`.
#' @aliases probability_distribution,ExMarkovParam-method
#'
#' @inheritParams probability_distribution
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#'
#' @examples
#' probability_distribution(CuadrasAugeExtMO2FParam(
#'   dim = 5L, lambda = 0.05, rho = 0.4), 0.3)
#' probability_distribution(AlphaStableExtMO2FParam(
#'   dim = 5L, lambda = 0.05, rho = 0.4), 0.3)
#' probability_distribution(PoissonExtMO2FParam(
#'   dim = 5L, lambda = 0.05, rho = 0.4), 0.3)
#' probability_distribution(ExponentialExtMO2FParam(
#'   dim = 5L, lambda = 0.05, rho = 0.4), 0.3)
#'
#' @section Probability distribution:
#' The probability of \eqn{j > i} portfolio items being defaulted at time
#' \eqn{t > s} conditioned on \eqn{i} portfolio items being defaulted at time
#' \eqn{s} is
#'
#' \deqn{
#'   \mathbb{P}(Z_t = j \mid Z_s = i)
#'     = \delta_{i}^\top \operatorname{e}^{(t-s) Q} \delta_{j} .
#' }
#'
#' @importFrom expm expm
#' @importFrom checkmate qassert
#' @importFrom purrr map reduce
#' @include utils.R
#'
#' @export
setMethod("probability_distribution", "ExMarkovParam",
  function(object, times, ...,
      method = c("default", "ExMarkovParam", "CalibrationParam")) {
    method <- match.arg(method)
    if (isTRUE("default" == method)) {
      method <- "ExMarkovParam"
    }
    if (!isTRUE("ExMarkovParam" == method)) {
      out <- callNextMethod(object, times, ..., method = method)
    } else {
      qassert(times, "N+[0,)")
      out <- map(times, ~expm(. * object@ex_qmatrix)[1, ]) %>%
        reduce(cbind) %>%
        `dimnames<-`(NULL) %>%
        matrix(nrow = getDimension(object) + 1L, ncol = length(times))
    }

    out
  })


#' @describeIn ExMarkovParam-class Display the object.
#' @aliases show,ExMarkovParam-method
#'
#' @export
setMethod("show", "ExMarkovParam",
  function(object) {
    cat("An object of class \"ExMarkovParam\"\n")
    cat(sprintf("Dimension: %i\n", getDimension(object)))
    cat("Generator matrix:\n")
    print.table(getExQMatrix(object), zero.print = "")

    invisible(NULL)
  })
