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
setGeneric("simulate_adcp",
  function(object, times, ...) {
    standardGeneric("simulate_adcp")
  })

#' @describeIn CalibrationParam-class
#'    simulates the default times \eqn{(\tau_1, \ldots, \tau_d)} and returns a
#'    matrix `x` with `nrow(x) == n_sim` and `ncol(x) == dim(object)` if
#'    `dim(object) > 1L` and a vector `x` with `length(x) == n_sim` otherwise.
#'
#' @param object [CalibrationParam-class]-object.
#' ... Pass-through parameters.
#'
#' @export
setGeneric("simulate_dt",
  function(object, ...) {
    standardGeneric("simulate_dt")
  })

#' @describeIn CalibrationParam-class
#'    simulates the default times \eqn{(\tau_1, \ldots, \tau_d)} and returns a
#'    matrix `x` with `nrow(x) == n_sim` and `ncol(x) == dim(object)` if
#'    `dim(object) > 1L` and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_adcp,CalibrationParam-method
#'
#' @examples
#' parm <- ExMarkovParam(
#'  ex_qmatrix = matrix(
#'    c(-0.07647059, 0, 0, 0.05294118, -0.05, 0, 0.02352941, 0.05, 0),
#'    nrow = 3L, ncol = 3L))
#' simulate_adcp(parm, 1, n_sim = 5e1)
#' simulate_adcp(parm, seq(0, 5, by = 0.25), n_sim = 5e1)
#'
#' @importFrom stats rexp
#' @export
setMethod("simulate_adcp", "CalibrationParam",
  function(object, times, ...) {
    tmp <- simulate_dt(object, ...)
    if (!is.matrix(tmp)) {
      tmp <- matrix(tmp, ncol = dim(object))
    }
    simplify2vector(dt2adcp(tmp, times))
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
    }

    simplify2vector(out)
  })

#' @describeIn ExMOParam-class
#'    simulates the default times \eqn{(\tau_1, \ldots, \tau_d)} and returns a
#'    matrix `x` with `nrow(x) == n_sim` and `ncol(x) == dim(object)` if
#'    `dim(object) > 1L` and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_dt,ExMOParam-method
#'
#' @inheritParams simulate_dt
#' @param method Simulation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param n_sim Number of samples.
#'
#' @examples
#' parm <- ExMOParam(ex_intensities = c(0.02647059, 0.02352941))
#' simulate_dt(parm, n_sim = 5e1)
#'
#' @importFrom rmo rexmo_markovian
#' @export
setMethod("simulate_dt", "ExMOParam",
  function(object, ...,
      method = c("default", "ExMOParam", "ExMarkovParam"), n_sim = 1e4) {
    method <- match.arg(method)
    if (isTRUE("default" == method)) {
      method <- "ExMOParam"
    }

    if (isTRUE("ExMOParam" == method)) {
      out <- rexmo_markovian(n_sim, object@dim, object@ex_intensities)
      out <- simplify2vector(out)
    } else {
      out <- callNextMethod(object, ..., n_sim = n_sim)
    }

    out
  })

#' @describeIn ExtGaussian2FParam-class
#'    simulates the default times \eqn{(\tau_1, \ldots, \tau_d)} and returns a
#'    matrix `x` with `nrow(x) == n_sim` and `ncol(x) == dim(object)` if
#'    `dim(object) > 1L` and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_dt,ExtGaussian2FParam-method
#'
#' @inheritParams simulate_dt
#' @param method Simulation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param n_sim Number of samples.
#'
#' @examples
#' parm <- ExtGaussian2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
#' simulate_dt(parm, n_sim = 5e1)
#'
#' @importFrom stats qexp
#' @importFrom copula normalCopula rCopula
setMethod("simulate_dt", "ExtGaussian2FParam",
  function(object, ...,
      method = c("default", "ExtGaussian2FParam"), n_sim = 1e4) {
    method <- match.arg(method)
    out <- qexp(
      rCopula(
        n_sim,
        normalCopula(param = object@nu, dim = object@dim, dispstr = "ex")
      ),
      rate = object@lambda, lower.tail = FALSE
    )

    simplify2vector(out)
  })

#' @describeIn ExtArch2FParam-class
#'    simulates the default times \eqn{(\tau_1, \ldots, \tau_d)} and returns a
#'    matrix `x` with `nrow(x) == n_sim` and `ncol(x) == dim(object)` if
#'    `dim(object) > 1L` and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_dt,ExtArch2FParam-method
#'
#' @inheritParams simulate_dt
#' @param method Simulation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param n_sim Number of samples.
#'
#' @examples
#' parm <- FrankExtArch2FParam(dim = 5L, lambda = 8e-2, rho = 4e-1)
#' simulate_dt(parm, n_sim = 5e1)
#'
#' @importFrom stats qexp
#' @importFrom copula rCopula
setMethod("simulate_dt", "ExtArch2FParam",
  function(object, ...,
      method = c("default", "ExtArch2fParam"), n_sim = 1e4) {
    method <- match.arg(method)
    out <- qexp(
      rCopula(n_sim, object@copula),
      rate = object@lambda, lower.tail = !object@survival
    )

    simplify2vector(out)
  })

#' @describeIn H2ExMarkovParam-class
#'    simulates the default times \eqn{(\tau_1, \ldots, \tau_d)} and returns a
#'    matrix `x` with `nrow(x) == n_sim` and `ncol(x) == dim(object)` if
#'    `dim(object) > 1L` and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_dt,H2ExMarkovParam-method
#'
#' @inheritParams simulate_dt
#'
#' @examples
#' parm <- ExponentialH2ExtMO3FParam(
#'   partition = list(1:2, 3:6, 7:8), fraction = 0.4,
#'   lambda = 8e-2, rho = c(0.2, 0.7))
#' simulate_dt(parm, n_sim = 5e1)
setMethod("simulate_dt", "H2ExMarkovParam",
  function(object, ...) {
    fraction <- object@fraction
    tmp0 <- simulate_dt(object@models[[1]], ...)
    tmp1 <- reduce(map(object@models[-1], simulate_dt, ...), cbind)
    out <- pmin(1 / fraction * tmp0, 1 / (1 - fraction) * tmp1)

    simplify2vector(out)
  })

#' @describeIn H2ExtGaussian3FParam-class
#'    simulates the default times \eqn{(\tau_1, \ldots, \tau_d)} and returns a
#'    matrix `x` with `nrow(x) == n_sim` and `ncol(x) == dim(object)` if
#'    `dim(object) > 1L` and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_dt,H2ExtGaussian3FParam-method
#'
#' @inheritParams simulate_dt
#' @param n_sim Number of samples.
#'
#' @examples
#' parm <- H2ExtGaussian3FParam(
#'   partition = list(1:2, 3:6, 7:8),
#'   lambda = 8e-2, rho = c(0.2, 0.7))
#' simulate_dt(parm, n_sim = 5e1)
#'
#' @importFrom purrr walk
#' @importFrom copula normalCopula rCopula P2p p2P
#' @importFrom stats qexp
setMethod("simulate_dt", "H2ExtGaussian3FParam",
  function(object, ..., n_sim = 1e4) {
    d <- getDimension(object)
    nu <- getNu(object)
    corr <- p2P(nu[[1]], d = d)
    walk(object@partition, ~{
      corr[.x, .x] <<- p2P(nu[[2]], d = length(.x))
    })

    out <- qexp(
      rCopula(
        n_sim,
        normalCopula(param = P2p(corr), dim = d, dispstr = "un")
      ),
      rate = getLambda(object), lower.tail = FALSE
    )

    simplify2vector(out)
  })

#' @describeIn H2ExtArch3FParam-class
#'    simulates the default times \eqn{(\tau_1, \ldots, \tau_d)} and returns a
#'    matrix `x` with `nrow(x) == n_sim` and `ncol(x) == dim(object)` if
#'    `dim(object) > 1L` and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_dt,H2ExtArch3FParam-method
#'
#' @inheritParams simulate_dt
#' @param n_sim Number of samples.
#'
#' @examples
#' parm <- FrankH2ExtArch3FParam(
#'   partition = list(1:2, 3:6, 7:8),
#'   lambda = 8e-2, rho = c(0.2, 0.7))
#' simulate_dt(parm, n_sim = 5e1)
#'
#' @importFrom stats qexp
#' @importFrom copula rCopula
setMethod("simulate_dt", "H2ExtArch3FParam",
  function(object, ..., n_sim = 1e4) {
    out <- qexp(
      rCopula(n_sim, object@copula),
      rate = object@lambda, lower.tail = !object@survival
    )

    simplify2vector(out)
  })
