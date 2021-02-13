#' @include allClass.R allGeneric.R allProbability.R
NULL

#' Expected values for the calibration parameter
#'
#' Calculates expected values for the average default counting process \eqn{L}.
#'
#' @inheritParams probability_distribution
#' @param g Transformation function
#' @param ... Further arguments to `g`
#'
#' @details
#' Calculates for a function \eqn{g} and the average default counting process
#' \eqn{L} the expectation
#' \deqn{
#'   \mathbb{E}[g(L_t)] .
#' }
#' For a portfolio CDS choose \eqn{g(x) = (1 - R) x} and for a CDO tranche with
#' attachment points \eqn{l < u} and choose
#' \eqn{g(x) = min{\{ \max{\{ (1 - R) x - l, 0 \}}, u - l \}}}, where \eqn{R} is
#' the recovery rate.
#'
#' @docType methods
#' @export
setGeneric("expected_value",
  function(object, times, g, ...) {
    standardGeneric("expected_value")
  })

#' @rdname expected_value
#'
#' @param recovery_rate The recovery rate of the portfolio CDS/CDO
#' @param ... Further arguments
#'
#' @export
setGeneric("expected_pcds_loss",
  function(object, times, recovery_rate, ...) {
    standardGeneric("expected_pcds_loss")
  })

#' @rdname expected_value
#'
#' @param lower Lower attachment point of the CDO tranche
#' @param upper Upper attachment point of the CDO tranche
#'
#' @export
setGeneric("expected_cdo_loss",
  function(object, times, recovery_rate, lower, upper, ...) {
    standardGeneric("expected_cdo_loss")
  })

#' @rdname expected_value
#' @aliases expected_value,CalibrationParam
#'
#' @examples
#' expected_value(CuadrasAugeExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3,
#'   function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2))
#' expected_value(AlphaStableExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3,
#'   function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2))
#' expected_value(PoissonExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3,
#'   function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2))
#' expected_value(ExponentialExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3,
#'   function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2))
#'
#' @importFrom checkmate qassert assert_function
#' @export
setMethod("expected_value", "CalibrationParam",
  function(object, times, g, ...) {
    qassert(times, "N+[0,)")
    assert_function(g)
    mu <- sapply(0:object@dim, function(k) g(k / object@dim, ...))

    as.vector(t(probability_distribution(object, times)) %*% mu)
  })

#' @rdname expected_value
#' @aliases expected_pcds_loss,CalibrationParam
#'
#' @examples
#' expected_pcds_loss(CuadrasAugeExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4)
#' expected_pcds_loss(AlphaStableExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4)
#' expected_pcds_loss(PoissonExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4)
#' expected_pcds_loss(ExponentialExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4)
#'
#' @importFrom checkmate qassert
#' @export
setMethod("expected_pcds_loss", "CalibrationParam",
  function(object, times, recovery_rate, ...) {
    qassert(recovery_rate, "N1[0,1]")
    expected_value(object, times,
      function(k) {
        (1 - recovery_rate) * k
      })
  })

#' @rdname expected_value
#'
#' @param method Choice of method (if available)
#'
#' @examples
#' expected_pcds_loss(CuadrasAugeExtMO2FParam(dim = 75, lambda = 0.05, rho = 0.4),
#'   times = 0.25, recovery_rate = 0.4)
#' expected_pcds_loss(CuadrasAugeExtMO2FParam(dim = 75, lambda = 0.05, rho = 0.4),
#'   times = 0.25, recovery_rate = 0.4, method = "fallback")
#'
#' @importFrom stats pexp
#' @importFrom checkmate qassert
#' @export
setMethod("expected_pcds_loss", "ExtMO2FParam",
  function(object, times, recovery_rate, method = c("default", "fallback"), ...) {
    method <- match.arg(method)
    qassert(times, "N+[0,)")
    qassert(recovery_rate, "N1[0,1]")
    if (!isTRUE("default" == method)) {
      return(callNextMethod(object, times, recovery_rate, ...))
    }

    (1 - recovery_rate) * pexp(times, rate = object@lambda)
  })

#' @rdname expected_value
#' @aliases expected_pcds_loss,ExtGaussian2FParam
#'
#' @examples
#' expected_pcds_loss(ExtGaussian2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   times = 0.25, recovery_rate = 0.4)
#' expected_pcds_loss(ExtGaussian2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   times = 0.25, recovery_rate = 0.4, method = "fallback")
#'
#' @importFrom stats pexp
#' @importFrom checkmate qassert
#' @export
setMethod("expected_pcds_loss", "ExtGaussian2FParam",
  function(object, times, recovery_rate, method = c("default", "fallback"), ...) {
    method <- match.arg(method)
    qassert(times, "N+[0,)")
    qassert(recovery_rate, "N1[0,1]")
    if (!isTRUE("default" == method)) {
      return(callNextMethod(object, times, recovery_rate, ...))
    }

    (1 - recovery_rate) * pexp(times, rate = object@lambda)
  })

#' @rdname expected_value
#' @aliases expected_pcds_loss,FrankExtArch2FParam
#'
#' @examples
#' expected_pcds_loss(FrankExtArch2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   times = 0.25, recovery_rate = 0.4)
#' expected_pcds_loss(FrankExtArch2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   times = 0.25, recovery_rate = 0.4, method = "fallback")
#'
#' @importFrom stats pexp
#' @importFrom checkmate qassert
#' @export
setMethod("expected_pcds_loss", "FrankExtArch2FParam",
  function(object, times, recovery_rate, method = c("default", "fallback"), ...) {
    method <- match.arg(method)
    qassert(times, "N+[0,)")
    qassert(recovery_rate, "N1[0,1]")
    if (!isTRUE("default" == method)) {
      return(callNextMethod(object, times, recovery_rate, ...))
    }

    (1 - recovery_rate) * pexp(times, rate = object@lambda)
  })

#' @rdname expected_value
#' @aliases expected_cdo_loss,CalibrationParam
#'
#' @examples
#' expected_cdo_loss(CuadrasAugeExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4, 0.1, 0.2)
#' expected_cdo_loss(AlphaStableExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4, 0.1, 0.2)
#' expected_cdo_loss(PoissonExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4, 0.1, 0.2)
#' expected_cdo_loss(ExponentialExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3, 0.4, 0.1, 0.2)
#'
#' @importFrom checkmate qassert assert_numeric
#' @export
setMethod("expected_cdo_loss", "CalibrationParam",
  function(object, times, recovery_rate, lower, upper, ...) {
    qassert(recovery_rate, "N1[0,1]")
    qassert(lower, "N1[0,1]")
    assert_numeric(upper, lower = lower, upper = 1)
    expected_value(object, times,
      function(k) {
        pmin(pmax((1 - recovery_rate) * k - lower, 0), upper - lower)
      })
  })

#' @rdname expected_value
#' @aliases expected_cdo_loss,ExtGaussian2FParam
#'
#' @examples
#' expected_cdo_loss(ExtGaussian2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   times = 0.25, recovery_rate = 0.4, lower = 0.1, upper = 0.2)
#' expected_cdo_loss(ExtGaussian2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   times = 0.25, recovery_rate = 0.4, lower = 0.1, upper = 0.2, method = "fallback")
#'
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats pexp qnorm
#' @importFrom checkmate qassert assert_numeric
#' @export
setMethod("expected_cdo_loss", "ExtGaussian2FParam",
  function(
      object, times, recovery_rate, lower, upper,
      method = c("default", "fallback"), ...) {
    method <- match.arg(method)
    qassert(times, "N+[0,)")
    qassert(recovery_rate, "N1[0,1]")
    qassert(lower, "N1[0,1]")
    assert_numeric(upper, lower = lower, upper = 1)
    if (!isTRUE("default" == method)) {
      return(callNextMethod(object, times, recovery_rate, lower, upper, ...))
    }

    corr <- matrix(c(1, rep(-sqrt(1 - object@nu), 2), 1), nrow=2, ncol=2)
    times <- qnorm(pexp(times, rate = object@lambda))
    sapply(times, function(t) {
      left <- pnorm(t)
      if (lower > 0) {
        left <- pmvnorm(
          lower = rep(-Inf, 2),
          upper = c(
            -qnorm(pmin(lower / (1 - recovery_rate), 1)),
            t
          ),
          corr = corr
        )
      }
      right <- pmvnorm(
        lower = rep(-Inf, 2),
        upper = c(
          -qnorm(pmin(upper / (1 - recovery_rate), 1)),
          t
        ),
        corr = corr
      )

      (1 - recovery_rate) * as.numeric(left - right)
    })
  })
