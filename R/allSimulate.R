#' @include allClass.R allGeneric.R
NULL

#' Simulate the default counting process
#'
#' Simulation of the average default counting process \eqn{L} for a given
#' time grid.
#'
#' @param object Calibration parameter object
#' @param n Number of samples
#' @param times The times for which the process is supposed to be simulated
#' @param ... Further parameters
#'
#' @export
setGeneric("simulate_param",
  function(object, n, times, ...) {
    standardGeneric("simulate_param")
  })

#' @rdname simulate_param
#' @aliases simulate_param,ExMarkovParam
#'
#' @examples
#' simulate_param(ExMarkovParam(), 1e1, seq(0, 5, by = 0.25))
#'
#' @importFrom stats rexp
#' @export
setMethod("simulate_param", "ExMarkovParam",
  function(object, n, times, ...) {
    out <- matrix(nrow = n, ncol = object@dim)
    for (k in 1:n) {
      state <- 0
      time <- 0
      while (state != object@dim) {
        wt <- rexp(1, rate = -object@qmatrix[1+state, 1+state])
        out[k, (1+state):object@dim]
        state <- state + sample.int(n = object@dim-state, size = 1, replace = FALSE,
                        prob = object@qmatrix[1+state, (2+state):(object@dim+1)])
      }
    }

    apply(out, 1, function(x) sapply(times, function(t) mean(x <= t)))
  })

#' @rdname simulate_param
#' @aliases simulate_param,ExMOParam
#'
#' @examples
#' simulate_param(ExMOParam(), 1e1, seq(0, 5, by = 0.25))
#'
#' @importFrom rmo rexmo_markovian
#' @export
setMethod("simulate_param", "ExMOParam",
  function(object, n, times, ...) {
    out <- rexmo_markovian(n, object@dim, object@ex_intensities)

    apply(out, 1, function(x) sapply(times, function(t) mean(x <= t)))
  })

#' @rdname simulate_param
#' @aliases simulate_param,ExtGaussian2FParam
#'
#' @examples
#' simulate_param(ExtGaussian2FParam(dim = 5), 1e1, seq(0, 5, by = 0.25))
#'
#' @importFrom stats qexp
#' @importFrom copula normalCopula rCopula
setMethod("simulate_param", "ExtGaussian2FParam",
  function(object, n, times, ...) {
    out <- qexp(
      rCopula(
        n,
        normalCopula(param = object@nu, dim = object@dim, dispstr = "ex")
      ),
      rate = object@lambda, lower.tail = FALSE
    )

    apply(out, 1, function(x) sapply(times, function(t) mean(x <= t)))
  })

#' @rdname simulate_param
#' @aliases simulate_param,FrankExtArch2FParam
#'
#' @examples
#' simulate_param(FrankExtArch2FParam(dim = 5), 1e1, seq(0, 5, by = 0.25))
#'
#' @importFrom stats qexp
#' @importFrom copula frankCopula rCopula
setMethod("simulate_param", "FrankExtArch2FParam",
  function(object, n, times, ...) {
    out <- qexp(
      rCopula(
        n,
        frankCopula(param = object@nu, dim = object@dim)
      ),
      rate = object@lambda, lower.tail = FALSE
    )

    apply(out, 1, function(x) sapply(times, function(t) mean(x <= t)))
  })
