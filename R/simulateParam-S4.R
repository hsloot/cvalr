#' @include allClass-S4.R allGeneric-S4.R utils.R
NULL

#' @describeIn CalibrationParam-class
#'   simulates the default counting process \eqn{L} and returns a matrix `x` with
#'   `nrow(x) == n_sim` and `ncol(x) == length(times)` if `length(times) > 1L`
#'   and a vector `x` with `length(x) == n_sim` otherwise.
#'
#' @param object Calibration parameter object.
#' @param times The times for which the process is supposed to be simulated.
#' @param ... Pass-through parameter.
#'
#' @export
setGeneric("simulate_param",
  function(object, times, ...) {
    standardGeneric("simulate_param")
  })

#' @describeIn ExMarkovParam-class
#'   simulates the default counting process \eqn{L} and returns a matrix `x` with
#'   `nrow(x) == n_sim` and `ncol(x) == length(times)` if `length(times) > 1L`
#'   and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_param,ExMarkovParam-method
#'
#' @inheritParams simulate_param
#' @param n_sim Number of samples.
#'
#' @examples
#' simulate_param(ExMarkovParam(), 1e1, seq(0, 5, by = 0.25))
#'
#' @importFrom stats rexp
#' @export
setMethod("simulate_param", "ExMarkovParam",
  function(object, times, ..., n_sim = 1e4) {
    tmp <- matrix(nrow = n_sim, ncol = object@dim)
    for (k in 1:n_sim) {
      state <- 0
      time <- 0
      while (state != object@dim) {
        wt <- rexp(1, rate = -object@qmatrix[1+state, 1+state])
        time <- time + wt
        tmp[k, (1+state):object@dim] <- time
        state <- state + sample.int(n = object@dim-state, size = 1, replace = FALSE,
                        prob = object@qmatrix[1+state, (2+state):(object@dim+1)])
      }
    }
    out <- t(apply(tmp, 1, function(x) sapply(times, function(t) mean(x <= t))))

    simplify2vector(out)
  })

#' @describeIn ExMOParam-class
#'   simulates the default counting process \eqn{L} and returns a matrix `x` with
#'   `nrow(x) == n_sim` and `ncol(x) == length(times)` if `length(times) > 1L`
#'   and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_param,ExMOParam-method
#'
#' @inheritParams simulate_param
#' @param method Simulation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param n_sim Number of samples.
#'
#' @examples
#' simulate_param(ExMOParam(), 1e1, seq(0, 5, by = 0.25))
#'
#' @importFrom rmo rexmo_markovian
#' @export
setMethod("simulate_param", "ExMOParam",
  function(object, times, ...,
      method = c("default", "ExMOParam", "ExMarkovParam"), n_sim = 1e4) {
    method <- match.arg(method)
    if (isTRUE("default" == method)) {
      method <- "ExMOParam"
    }

    if (isTRUE("ExMOParam" == method)) {
      tmp <- rexmo_markovian(n_sim, object@dim, object@ex_intensities)
      out <- t(apply(tmp, 1, function(x) sapply(times, function(t) mean(x <= t))))
    } else {
      out <- callNextMethod(object, times, ..., n_sim = n_sim)
    }

    simplify2vector(out)
  })

#' @describeIn ExtGaussian2FParam-class
#'   simulates the default counting process \eqn{L} and returns a matrix `x` with
#'   `nrow(x) == n_sim` and `ncol(x) == length(times)` if `length(times) > 1L`
#'   and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_param,ExtGaussian2FParam-method
#'
#' @inheritParams simulate_param
#' @param n_sim Number of samples.
#'
#' @examples
#' simulate_param(ExtGaussian2FParam(dim = 5), 1e1, seq(0, 5, by = 0.25))
#'
#' @importFrom stats qexp
#' @importFrom copula normalCopula rCopula
setMethod("simulate_param", "ExtGaussian2FParam",
  function(object, times, ..., n_sim = 1e4) {
    tmp <- qexp(
      rCopula(
        n_sim,
        normalCopula(param = object@nu, dim = object@dim, dispstr = "ex")
      ),
      rate = object@lambda, lower.tail = FALSE
    )
    out <- t(apply(tmp, 1, function(x) sapply(times, function(t) mean(x <= t))))

    simplify2vector(out)
  })

#' @describeIn ExtArch2FParam-class
#'   simulates the default counting process \eqn{L} and returns a matrix `x` with
#'   `nrow(x) == n_sim` and `ncol(x) == length(times)` if `length(times) > 1L`
#'   and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_param,FrankExtArch2FParam-method
#'
#' @inheritParams simulate_param
#' @param n_sim Number of samples.
#'
#' @examples
#' simulate_param(FrankExtArch2FParam(dim = 5), 1e1, seq(0, 5, by = 0.25))
#'
#' @importFrom stats qexp
#' @importFrom copula frankCopula rCopula
setMethod("simulate_param", "FrankExtArch2FParam",
  function(object, times, ..., n_sim = 1e4) {
    tmp <- qexp(
      rCopula(
        n_sim,
        frankCopula(param = object@nu, dim = object@dim)
      ),
      rate = object@lambda, lower.tail = FALSE
    )
    out <- t(apply(tmp, 1, function(x) sapply(times, function(t) mean(x <= t))))

    simplify2vector(out)
  })
