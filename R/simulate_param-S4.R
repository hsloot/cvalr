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
#' @param method Simulation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param n_sim Number of samples.
#'
#' @examples
#' parm <- ExMarkovParam(
#'  ex_qmatrix = matrix(
#'    c(-0.07647059, 0, 0, 0.05294118, -0.05, 0, 0.02352941, 0.05, 0),
#'    nrow = 3L, ncol = 3L))
#' simulate_param(parm, 1e1, seq(0, 5, by = 0.25), n_sim = 5e1)
#'
#' @importFrom stats rexp
#' @export
setMethod("simulate_param", "ExMarkovParam",
  function(object, times, ...,
      method = c("default", "ExMarkovParam"), n_sim = 1e4) {
    method <- match.arg(method)
    tmp <- matrix(nrow = n_sim, ncol = object@dim)
    for (k in 1:n_sim) {
      state <- 0
      time <- 0
      while (state != object@dim) {
        wt <- rexp(1, rate = -object@ex_qmatrix[1+state, 1+state])
        time <- time + wt
        tmp[k, (1+state):object@dim] <- time
        state <- state + sample.int(n = object@dim-state, size = 1, replace = FALSE,
                        prob = object@ex_qmatrix[1+state, (2+state):(object@dim+1)])
      }
    }
    out <- dt2adcp(tmp, times)

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
#' parm <- ExMOParam(ex_intensities = c(0.02647059, 0.02352941))
#' simulate_param(parm, 1e1, seq(0, 5, by = 0.25), n_sim = 5e1)
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
      out <- dt2adcp(tmp, times)
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
#' @param method Simulation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param n_sim Number of samples.
#'
#' @examples
#' parm <- ExtGaussian2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
#' simulate_param(parm, 1e1, seq(0, 5, by = 0.25), n_sim = 5e1)
#'
#' @importFrom stats qexp
#' @importFrom copula normalCopula rCopula
setMethod("simulate_param", "ExtGaussian2FParam",
  function(object, times, ...,
      method = c("default", "ExtGaussian2FParam"), n_sim = 1e4) {
    method <- match.arg(method)
    tmp <- qexp(
      rCopula(
        n_sim,
        normalCopula(param = object@nu, dim = object@dim, dispstr = "ex")
      ),
      rate = object@lambda, lower.tail = FALSE
    )
    out <- dt2adcp(tmp, times)

    simplify2vector(out)
  })

#' @describeIn ExtArch2FParam-class
#'   simulates the default counting process \eqn{L} and returns a matrix `x` with
#'   `nrow(x) == n_sim` and `ncol(x) == length(times)` if `length(times) > 1L`
#'   and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_param,ExtArch2FParam-method
#'
#' @inheritParams simulate_param
#' @param method Simulation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param n_sim Number of samples.
#'
#' @examples
#' parm <- FrankExtArch2FParam(dim = 5L, lambda = 8e-2, rho = 4e-1)
#' simulate_param(parm, 1e1, seq(0, 5, by = 0.25), n_sim = 5e1)
#'
#' @importFrom stats qexp
#' @importFrom copula rCopula
setMethod("simulate_param", "ExtArch2FParam",
  function(object, times, ...,
      method = c("default", "ExtArch2fParam"), n_sim = 1e4) {
    method <- match.arg(method)
    tmp <- qexp(
      rCopula(n_sim, object@copula),
      rate = object@lambda, lower.tail = !object@survival
    )
    out <- dt2adcp(tmp, times)

    simplify2vector(out)
  })
