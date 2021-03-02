#' @include allClass-S4.R allGeneric-S4.R simulate_param-S4.R utils.R
NULL

#' @describeIn CalibrationParam-class
#'   returns the probability vector for the average default count process \eqn{L}.
#'
#' @param object The calibration parameter object.
#' @param times Point-in-time.
#' @param ... Pass-through parameter.
#'
#' @section Probability distribution:
#' The probability vector of the *average default counting process* \eqn{L}
#' for certain times can be calculated with [probability_distribution()];
#' i.e. the values
#' \deqn{
#'   \mathbb{P}(L_t = k/d) , \quad k \in {\{ 0, \ldots, d \}} , \quad t \geq 0 .
#' }
#'
#' @export
setGeneric("probability_distribution",
  function(object, times, ...) {
    standardGeneric("probability_distribution")
  })

#' @describeIn ExMarkovParam-class
#'   returns the probability vector for the average default count process \eqn{L}.
#' @aliases probability_distribution,CalibrationParam-method
#'
#' @inheritParams probability_distribution
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param seed Numeric number (if not NULL, is used to set the seed prior to
#'   Monte-Carlo estimation of probability distribution).
#' @param sim_args List with pass-through parameters for [simulate_param()].
#'
#' @importFrom checkmate qassert assert_number assert_list
#' @export
setMethod("probability_distribution", "CalibrationParam",
  function(object, times, ...,
      method = c("default", "CalibrationParam"), seed = NULL, sim_args = NULL) {
    method <- match.arg(method)
    qassert(times, "N+[0,)")
    assert_number(seed, lower = 0, finite = TRUE, null.ok = TRUE)
    assert_list(sim_args, null.ok = TRUE)
    if (!is.null(seed)) {
      set.seed(seed)
    }
    x <- do.call(simulate_param,
           args = c(list(object = object, times = times), sim_args))
    if (!is.matrix(x)) {
      x <- as.matrix(x, ncol = 1L)
    }
    out <- adcp2epd(x, object@dim)

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
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#' probability_distribution(AlphaStableExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#' probability_distribution(PoissonExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#' probability_distribution(ExponentialExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#'
#' @importFrom expm expm
#' @importFrom checkmate qassert
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

#' @describeIn ExtGaussian2FParam-class
#'   returns the probability vector for the average default count process \eqn{L}.
#' @aliases probability_distribution,ExtGaussian2FParam-method
#'
#' @inheritParams probability_distribution
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#'
#' @examples
#' probability_distribution(ExtGaussian2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#'
#' @importFrom stats integrate pexp pnorm dnorm qnorm
#' @importFrom checkmate qassert
#' @export
setMethod("probability_distribution", "ExtGaussian2FParam",
  function(object, times, ...,
      method = c("default", "ExtGaussian2FParam", "CalibrationParam")) {
    method <- match.arg(method)
    if (isTRUE("default" == method)) {
      method <- "ExtGaussian2FParam"
    }
    if (!isTRUE("ExtGaussian2FParam" == method)) {
      out <- callNextMethod(object, times, ...)
    } else {
      qassert(times, "N+[0,)")
      times <- qnorm(pexp(times, rate = object@lambda))
      out <- outer(0:object@dim, times,
        Vectorize(function(k, t) {
          if (-Inf == t && 0 == k) {
            return(1)
          } else if (-Inf == t && 0 < k) {
            return(0)
          }
          int_res <- integrate(
            function(x) {
              ldp <- pnorm(
                (t - sqrt(object@nu) * x) / (sqrt(1 - object@nu)),
                log.p = TRUE, lower.tail = TRUE
              )
              lsp <- pnorm(
                (t - sqrt(object@nu) * x) / (sqrt(1 - object@nu)),
                log.p = TRUE, lower.tail = FALSE
              )
              v_multiply_binomial_coefficient(
                exp(k * ldp + (object@dim-k) * lsp) * dnorm(x), object@dim, k)
            }, lower = -Inf, upper = Inf, rel.tol = .Machine$double.eps^0.5
          )

          int_res$value
        }))
    }

    simplify2vector(out)
  })
