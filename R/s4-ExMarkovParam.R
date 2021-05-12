#' @include s4-CalibrationParam.R checkmate.R
NULL

#' Exchangeable Markovian calibration parameter
#'
#' Calibration parameter class for the general exchangeable model with a
#' Markovian *default counting process*.
#'
#' @slot qmatrix The \eqn{(d+1) \times (d+1)} Markov generator matrix of the
#'   default counting process.
#'
#' @details
#' The probability of \eqn{j > i} portfolio items being defaulted at time
#' \eqn{t > s} conditioned on \eqn{i} portfolio items being defaulted at time
#' \eqn{s} is
#'
#' \deqn{
#'   \mathbb{P}(Z_t = j \mid Z_s = i)
#'     = \delta_{i}^\top \operatorname{e}^{(t-s) Q} \delta_{j} .
#' }
#'
#' @export ExMarkovParam
ExMarkovParam <- setClass("ExMarkovParam", # nolint
  contains = "CalibrationParam",
  slots = c(ex_qmatrix = "matrix"))


setGeneric("getExQMatrix",
  function(object) {
    standardGeneric("getExQMatrix")
  })
setGeneric("setExQMatrix<-",
  function(object, value) {
    standardGeneric("setExQMatrix<-")
  })


setMethod("getExQMatrix", "ExMarkovParam",
  function(object) {
    object@ex_qmatrix
  })

setReplaceMethod("setExQMatrix", "ExMarkovParam",
  function(object, value) {
    assert_exqmatrix(value, min.rows = 1L, min.cols = 1L)

    dim <- nrow(value)-1
    setDimension(object) <- dim
    object@ex_qmatrix <- value

    invisible(object)
  })


setValidity("ExMarkovParam",
  function(object) {
    assert_exqmatrix(object@ex_qmatrix,
      nrows = object@dim+1L, ncols = object@dim + 1L)

    invisible(TRUE)
  })


#' @describeIn ExMarkovParam-class Constructor
#' @aliases initialize,ExMarkovParam-method
#' @aliases initialize,ExMarkovParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ex_qmatrix (Exchangeable) Q-matrix of the default counting process.
#'
#' @examples
#' ExMarkovParam(
#'  ex_qmatrix = matrix(
#'    c(-0.07647059, 0, 0, 0.05294118, -0.05, 0, 0.02352941, 0.05, 0),
#'    nrow = 3L, ncol = 3L))
setMethod("initialize", "ExMarkovParam",
  function(.Object, ex_qmatrix) { # nolint
    if (!missing(ex_qmatrix)) {
      setExQMatrix(.Object) <- ex_qmatrix
      validObject(.Object)
    }

    invisible(.Object)
  })


#' @describeIn ExMarkovParam-class
#'    simulates the default times \eqn{(\tau_1, \ldots, \tau_d)} and returns a
#'    matrix `x` with `nrow(x) == n_sim` and `ncol(x) == dim(object)` if
#'    `dim(object) > 1L` and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_dt,ExMarkovParam-method
#'
#' @inheritParams simulate_dt
#' @param method Simulation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param n_sim Number of samples.
#'
#' @examples
#' parm <- ExMarkovParam(
#'  ex_qmatrix = matrix(
#'    c(-0.07647059, 0, 0, 0.05294118, -0.05, 0, 0.02352941, 0.05, 0),
#'    nrow = 3L, ncol = 3L))
#' simulate_dt(parm, n_sim = 5e1)
#'
#' @importFrom stats rexp
#' @include utils.R
#' @export
setMethod("simulate_dt", "ExMarkovParam",
  function(object, ...,
      method = c("default", "ExMarkovParam"), n_sim = 1e4) {
    method <- match.arg(method)
    out <- matrix(nrow = n_sim, ncol = object@dim)
    for (k in 1:n_sim) {
      state <- 0
      time <- 0
      while (state != object@dim) {
        wt <- rexp(1, rate = -object@ex_qmatrix[1+state, 1+state])
        time <- time + wt
        out[k, (1+state):object@dim] <- time
        state <- state + sample.int(n = object@dim-state, size = 1, replace = FALSE,
                        prob = object@ex_qmatrix[1+state, (2+state):(object@dim+1)])
      }
      perm <- sample.int(n = object@dim, size = object@dim, replace = FALSE)
      out[k, perm] <- out[k, perm]
    }

    simplify2vector(out)
  })

#' @describeIn ExMarkovParam-class
#'   returns the probability vector for the average default count process \eqn{L}.
#' @aliases probability_distribution,ExMarkovParam-method
#'
#' @inheritParams probability_distribution
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#'
#' @examples
#' probability_distribution(CuadrasAugeExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.3)
#' probability_distribution(AlphaStableExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.3)
#' probability_distribution(PoissonExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.3)
#' probability_distribution(ExponentialExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.3)
#'
#' @importFrom expm expm
#' @importFrom checkmate qassert
#' @include utils.R
#' @export
setMethod("probability_distribution", "ExMarkovParam",
  function(object, times, ...,
      method = c("default", "ExMarkovParam", "CalibrationParam")) {
    method <- match.arg(method)
    if (isTRUE("default" == method)) {
      method <- "ExMarkovParam"
    }
    if (!isTRUE("ExMarkovParam" == method)) {
      out <- callNextMethod(object, times, ...)
    } else {
      qassert(times, "N+[0,)")
      out <- sapply(times, function(t) expm(t * object@ex_qmatrix)[1, ])
    }

    simplify2vector(out)
  })
