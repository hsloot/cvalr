#' @include s4-CalibrationParam.R checkmate.R
NULL

# nolint start
ERR_MSG_EXQMATRIX <- "`ex_qmatrix` must be upper-triangular Markov generator-matrix"
# nolint end

#' Exchangeable Markovian calibration parameter
#'
#' @description
#' [CalibrationParam-class]-class for the exchangeable Markovian *(average)
#' default counting process* model.
#'
#' @slot ex_qmatrix The \eqn{(d+1) \times (d+1)} Markov generator matrix (see
#'   [rmo::exQMatrix()]).
#'
#' @details
#' The model is defined by the assumption that the *(average) default counting
#' process* is Markovian on the state space \eqn{\{ 0, 1/d, \ldots, (d-1)/d, 1 \}}
#' with the generator matrix provided to the constructor.
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
    assert_exqmatrix(value, min.rows = 3L, min.cols = 3L)

    setDimension(object) <- nrow(value) - 1L
    object@ex_qmatrix <- value

    invisible(object)
  })


#' @importFrom checkmate test_matrix
#' @include RcppExports.R
setValidity("ExMarkovParam",
  function(object) {
    if (!(test_matrix(
          object@ex_qmatrix, mode = "numeric", any.missing = FALSE, all.missing = FALSE) &&
        is_exqmatrix(object@ex_qmatrix, tol = .Machine$double.eps^0.5))) {
      return(ERR_MSG_EXQMATRIX)
    }

    invisible(TRUE)
  })


#' @describeIn ExMarkovParam-class Constructor
#' @aliases initialize,ExMarkovParam-method
#' @aliases initialize,ExMarkovParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ex_qmatrix The \eqn{(d+1) \times (d+1)} Markov generator matrix (see
#'   [rmo::exQMatrix()]).
#'
#' @examples
#' ExMarkovParam()
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
#' @section Simulation:
#' The default times are sampled using the Markovian representation of the
#' *(average) default counting process*: The ordered version of the default-time
#' vector can be recorded while sampling the (average) default counting process.
#' This vector is uniformly shuffled to obtain a sample with the desired
#' distribution.
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
      method = c("default", "ExMarkovParam"), n_sim = 1e4L) {
    method <- match.arg(method)
    d <- getDimension(object)
    ex_qmatrix <- getExQMatrix(object)
    out <- matrix(nrow = n_sim, ncol = d)
    for (k in 1:n_sim) {
      state <- 0
      time <- 0
      while (state != d) {
        wt <- rexp(1, rate = -ex_qmatrix[1L+state, 1L+state])
        time <- time + wt
        out[k, (1L+state):d] <- time
        state <- state +
          sample.int(
            n = d-state, size = 1L, replace = FALSE,
            prob = ex_qmatrix[1L+state, (2L+state):(d+1L)])
      }
      perm <- sample.int(n = d, size = d, replace = FALSE)
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
    if (isTRUE("default" == method || "ExMarkovParam" == method)) {
      qassert(times, "N+[0,)")
      ex_qmatrix <- getExQMatrix(object)
      out <- map(times, ~expm(. * ex_qmatrix)[1L, ]) %>%
        reduce(cbind) %>%
        `dimnames<-`(NULL) %>%
        matrix(nrow = getDimension(object) + 1L, ncol = length(times))
    } else  {
      out <- callNextMethod(object, times, ..., method = method)
    }

    out
  })


#' @describeIn ExMarkovParam-class Display the object.
#' @aliases show,ExMarkovParam-method
#'
#' @export
setMethod("show", "ExMarkovParam",
  function(object) {
    cat(sprintf("An object of class %s\n", classLabel(class(object))))
    if (isTRUE(validObject(object, test = TRUE))) {
      cat(sprintf("Dimension: %i\n", getDimension(object)))
      cat("Generator matrix:\n")
      capture.output(print.table(getExQMatrix(object), zero.print = "")) %>%
        paste0("\t", .) %>%
        writeLines
    } else {
      cat("\t (invalid or not initialized)\n")
    }

    invisible(NULL)
  })
