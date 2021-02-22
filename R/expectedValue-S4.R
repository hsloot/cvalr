#' @include allClass-S4.R allGeneric-S4.R probabilityDistribution-S4.R
NULL

#' @describeIn CalibrationParam-class
#'   returns the expected value for the average default count process \eqn{L} at
#'   a specific time-point.
#'
#' @param g Transformation function.
#' @param ... Pass-through parameter.
#' @param pd_args Parameter for [probability_distribution()].
#'
#' @section Expected value:
#' The *expectated value* of the *average default counting process* \eqn{L}
#' under transformations can be calculated with [expected_value()]; i.e. the
#' value
#' \deqn{
#'   \mathbb{E}[g(L_t)] , \quad t \geq 0.
#' }
#' For a portfolio CDS choose \eqn{g(x) = (1 - R) x} and for a CDO tranche with
#' attachment points \eqn{l < u} and choose
#' \eqn{g(x) = min{\{ \max{\{ (1 - R) x - l, 0 \}}, u - l \}}}, where \eqn{R} is
#' the recovery rate.
#'
#' @export
setGeneric("expected_value",
  function(object, times, g, ..., pd_args = NULL) {
    standardGeneric("expected_value")
  })

#' @describeIn CalibrationParam-class
#'   returns the expected portfolio CDS loss for a specific time-point.
#'
#' @param recovery_rate The recovery rate of the portfolio CDS/CDO.
#'
#' @export
setGeneric("expected_pcds_loss",
  function(object, times, recovery_rate, ..., pd_args = NULL) {
    standardGeneric("expected_pcds_loss")
  })

#' @describeIn CalibrationParam-class
#'   returns the expected CDO loss for a specific time-point.
#'
#' @param lower Lower attachment point of the CDO tranche.
#' @param upper Upper attachment point of the CDO tranche.
#'
#' @export
setGeneric("expected_cdo_loss",
  function(object, times, recovery_rate, lower, upper, ..., pd_args = NULL) {
    standardGeneric("expected_cdo_loss")
  })

#' @describeIn CalibrationParam-class
#'   returns the expected portfolio CDS fair-value equation.
#'
#' @export
setGeneric("expected_pcds_equation",
  function(object, times, discount_factors, recovery_rate, coupon, upfront, ..., pd_args = NULL) {
    standardGeneric("expected_pcds_equation")
  })

#' @describeIn CalibrationParam-class
#'   returns the expected CDO fair-value equation.
#'
#' @export
setGeneric("expected_cdo_equation",
  function(object, times, discount_factors, recovery_rate, lower, upper, coupon, upfront, ..., pd_args = NULL) {
    standardGeneric("expected_cdo_equation")
  })

#' @rdname CalibrationParam-class
#' @aliases expected_value,CalibrationParam-method
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
  function(object, times, g, ..., pd_args = NULL) {
    qassert(times, "N+[0,)")
    assert_function(g)
    mu <- sapply(0:object@dim, function(k) g(k / object@dim, ...))

    as.vector(t(do.call(probability_distribution, args = c(list(object = object, times = times), pd_args))) %*% mu)
  })

#' @rdname CalibrationParam-class
#' @aliases expected_pcds_loss,CalibrationParam-method
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
  function(object, times, recovery_rate, ..., pd_args = NULL) {
    qassert(recovery_rate, "N1[0,1]")
    g <- function(k, recovery_rate) {
      (1 - recovery_rate) * k
    }
    expected_value(object, times, g, recovery_rate = recovery_rate, pd_args = pd_args)
  })

#' @describeIn ExtMO2FParam-class
#'   returns the expected portfolio CDS loss for a specific time-point.
#' @aliases expected_pcds_loss,ExtMO2FParam-method
#'
#' @inheritParams probability_distribution
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
  function(object, times, recovery_rate, method = c("default", "fallback"), ..., pd_args = NULL) {
    method <- match.arg(method)
    qassert(times, "N+[0,)")
    qassert(recovery_rate, "N1[0,1]")
    if (!isTRUE("default" == method)) {
      return(callNextMethod(object, times, recovery_rate, ...))
    }

    (1 - recovery_rate) * pexp(times, rate = object@lambda)
  })

#' @describeIn ExtGaussian2FParam-class
#'   returns the expected portfolio CDS loss for a specific time-point.
#' @aliases expected_pcds_loss,ExtGaussian2FParam-method
#'
#' @inheritParams probability_distribution
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
  function(object, times, recovery_rate, method = c("default", "fallback"), ..., pd_args = NULL) {
    method <- match.arg(method)
    qassert(times, "N+[0,)")
    qassert(recovery_rate, "N1[0,1]")
    if (!isTRUE("default" == method)) {
      return(callNextMethod(object, times, recovery_rate, ...))
    }

    (1 - recovery_rate) * pexp(times, rate = object@lambda)
  })

#' @describeIn ExtArch2FParam-class
#'   returns the expected portfolio CDS loss for a specific time-point.
#' @aliases expected_pcds_loss,FrankExtArch2FParam-method
#'
#' @inheritParams probability_distribution
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
  function(object, times, recovery_rate, method = c("default", "fallback"), ..., pd_args = NULL) {
    method <- match.arg(method)
    qassert(times, "N+[0,)")
    qassert(recovery_rate, "N1[0,1]")
    if (!isTRUE("default" == method)) {
      return(callNextMethod(object, times, recovery_rate, ...))
    }

    (1 - recovery_rate) * pexp(times, rate = object@lambda)
  })

#' @rdname CalibrationParam-class
#' @aliases expected_cdo_loss,CalibrationParam-method
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
  function(object, times, recovery_rate, lower, upper, ..., pd_args = NULL) {
    qassert(recovery_rate, "N1[0,1]")
    qassert(lower, "N1[0,1]")
    assert_numeric(upper, lower = lower, upper = 1)
    g <- function(k, recovery_rate, lower, upper) {
      pmin(pmax((1 - recovery_rate) * k - lower, 0), upper - lower)
    }
    expected_value(object, times, g,
      recovery_rate = recovery_rate, lower = lower, upper = upper,
      pd_args = pd_args)
  })

#' @describeIn ExtGaussian2FParam-class
#'   returns the expected CDO loss for a specific time-point.
#' @aliases expected_cdo_loss,ExtGaussian2FParam-method
#'
#' @inheritParams probability_distribution
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
      method = c("default", "fallback"), ..., pd_args = NULL) {
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

#' @rdname CalibrationParam-class
#' @aliases expected_pcds_equation,CalibrationParam-method
#'
#' @param discount_factors Discount factors for `times`
#' @param coupon Running coupon
#' @param upfront Upfront payment
#'
#' @examples
#' expected_pcds_equation(
#'   ExtGaussian2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   times = seq(0.25, 5, by = 0.25), discount_factors = rep(1, 20),
#'   recovery_rate = 0.4, coupon = 0.08, upfront = 0
#' )
#' expected_pcds_equation(
#'   ExtGaussian2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   times = seq(0.25, 5, by = 0.25), discount_factors = rep(1, 20),
#'   recovery_rate = rep(0.4, 4), coupon = rep(0.08, 4), upfront = rep(0, 4)
#' )
#'
#' @export
setMethod("expected_pcds_equation", "CalibrationParam",
function(object, times, discount_factors, recovery_rate, coupon, upfront, ..., pd_args = NULL) {
  qassert(coupon, "N+")
  qassert(upfront, "N+")
  expected_losses <- mapply(
    expected_pcds_loss,
    recovery_rate = recovery_rate,
    MoreArgs = c(list(object = object, times = times), list(...)),
    SIMPLIFY = FALSE
  )
  mapply(
    portfolio_cds_equation,
    expected_losses = expected_losses,
    recovery_rate = recovery_rate,
    coupon = coupon,
    upfront = upfront,
    MoreArgs = list(
      times = times,
      discount_factors = discount_factors
    ))
})

#' @rdname CalibrationParam-class
#' @aliases expected_cdo_equation,CalibrationParam-method
#'
#' @examples
#' expected_cdo_equation(
#'   ExtGaussian2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   times = seq(0.25, 5, by = 0.25), discount_factors = rep(1, 20),
#'   recovery_rate = 0.4, lower = c(0, 0.1, 0.2, 0.35),
#'   upper = c(0.1, 0.2, 0.35, 1), coupon = 0.08, upfront = 0
#' )
#'
#' @export
setMethod("expected_cdo_equation", "CalibrationParam",
function(object, times, discount_factors, recovery_rate, coupon, upfront, ..., pd_args = NULL) {
  qassert(coupon, "N+")
  qassert(upfront, "N+")
  expected_losses <- mapply(
    expected_cdo_loss,
    recovery_rate = recovery_rate,
    lower = lower, upper = upper,
    MoreArgs = c(list(object = object, times = times), list(...)),
    SIMPLIFY = FALSE
  )
  mapply(
    cdo_equation,
    expected_losses = expected_losses,
    lower = lower,
    upper = upper,
    coupon = coupon,
    upfront = upfront,
    MoreArgs = list(
      times = times,
      discount_factors = discount_factors
    ))
})
