#' @importFrom methods setClass setValidity setGeneric setMethod
#'   setReplaceMethod validObject is as new callNextMethod
#' @include checkmate.R
NULL

#' Virtual super-class for calibration parameters
#'
#' A virtual superclass for all calibration parameters for multivariate
#' (portfolio) credit models.
#'
#' @slot dim The dimension (number of portfolio items)
#'
#' @details
#' It is assumed that the expected value of the credit derivative
#' payment streams can be represented as a function of *expected losses* and
#' other non-model dependent quantities.
#' *Expected losses* are the expected values of the average default counting
#' process \eqn{L} under a transformation \eqn{g} on times \eqn{t_1, \ldots, t_n},
#' i.e. \eqn{E[g(L_{t_1})], \ldots, E[g(L_{t_n})]}.
#'
#' A `CalibrationParam` object is supposed to be used to abstract the calculation
#' of these *expected losses* in functions which calculate the fair value of the
#' option with the asset pricing theorem.
#'
#' @seealso [ExMarkovParam-class] [ExMOParam-class] [ExtMOParam-class]
#'   [CuadrasAugeExtMO2FParam-class] [AlphaStableExtMO2FParam-class]
#'   [PoissonExtMO2FParam-class] [ExponentialExtMO2FParam-class]
#'
#' @export
setClass("CalibrationParam", # nolint
  contains = "VIRTUAL",
  slots = c(dim = "integer"))


setGeneric("getDimension",
  function(object) {
    standardGeneric("getDimension")
  })
setGeneric("setDimension<-",
  function(object, value) {
    standardGeneric("setDimension<-")
  })


setMethod("getDimension", "CalibrationParam",
  function(object) {
    object@dim
  })

#' @importFrom checkmate qassert
setReplaceMethod("setDimension", "CalibrationParam",
  function(object, value) {
    qassert(value, "X1(0,)")
    object@dim <- as.integer(value)

    invisible(object)
  })


#' @importFrom checkmate qassert
setValidity("CalibrationParam",
  function(object) {
    qassert(object@dim, "I1(0,)")

    invisible(TRUE)
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
#' @include utils.R
#' @export
setMethod("simulate_adcp", "CalibrationParam",
  function(object, times, ...) {
    tmp <- simulate_dt(object, ...)
    if (!is.matrix(tmp)) {
      tmp <- matrix(tmp, ncol = dim(object))
    }
    simplify2vector(dt2adcp(tmp, times))
  })

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
#' @param sim_args List with pass-through parameters for [simulate_adcp()].
#'
#' @importFrom checkmate qassert assert_number assert_list
#' @include utils.R
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
    x <- do.call(simulate_adcp,
           args = c(list(object = object, times = times), sim_args))
    if (!is.matrix(x)) {
      x <- as.matrix(x, ncol = 1L)
    }
    out <- adcp2epd(x, object@dim)

    simplify2vector(out)
  })

#' @describeIn CalibrationParam-class
#'   returns the expected value for the average default count process \eqn{L} at
#'   a specific time-point.
#'
#' @param g Transformation function.
#' @param ... Pass-through parameter.
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
  function(object, times, g, ...) {
    standardGeneric("expected_value")
  })

#' @rdname CalibrationParam-class
#' @aliases expected_value,CalibrationParam-method
#'
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param pd_args Parameter for [probability_distribution()].
#'
#' @examples
#' expected_value(CuadrasAugeExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.3,
#'   function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2))
#' expected_value(AlphaStableExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.3,
#'   function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2))
#' expected_value(PoissonExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.3,
#'   function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2))
#' expected_value(ExponentialExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.3,
#'   function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2))
#'
#' @importFrom checkmate qassert assert_function assert_list
#' @export
setMethod("expected_value", "CalibrationParam",
  function(object, times, g, ...,
      method = c("default", "CalibrationParam"), pd_args = NULL) {
    method <- match.arg(method)
    qassert(times, "N+[0,)")
    assert_function(g)
    assert_list(pd_args, null.ok = TRUE)
    mu <- sapply(0:object@dim, function(k) g(k / object@dim, ...))

    probs <- do.call(probability_distribution,
                  args = c(list(object = object, times = times), pd_args))

    as.vector(t(probs) %*% mu)
  })

#' @describeIn CalibrationParam-class
#'   returns the expected portfolio CDS loss for a specific time-point.
#'
#' @param recovery_rate The recovery rate of the portfolio CDS/CDO.
#'
#' @export
setGeneric("expected_pcds_loss",
  function(object, times, recovery_rate, ...) {
    standardGeneric("expected_pcds_loss")
  })

#' @rdname CalibrationParam-class
#' @aliases expected_pcds_loss,CalibrationParam-method
#'
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#'
#' @examples
#' expected_pcds_loss(CuadrasAugeExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.3, 0.4)
#' expected_pcds_loss(AlphaStableExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.3, 0.4)
#' expected_pcds_loss(PoissonExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.3, 0.4)
#' expected_pcds_loss(ExponentialExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.3, 0.4)
#'
#' @importFrom checkmate qassert
#' @export
setMethod("expected_pcds_loss", "CalibrationParam",
  function(object, times, recovery_rate, ...,
      method = c("default", "CalibrationParam")) {
    method <- match.arg(method)
    qassert(recovery_rate, "N1[0,1]")
    g <- function(k, recovery_rate) {
      (1 - recovery_rate) * k
    }

    expected_value(object, times, g, recovery_rate = recovery_rate, ...)
  })

#' @describeIn CalibrationParam-class
#'   returns the expected CDO loss for a specific time-point.
#'
#' @param lower Lower attachment point of the CDO tranche.
#' @param upper Upper attachment point of the CDO tranche.
#'
#' @export
setGeneric("expected_cdo_loss",
  function(object, times, recovery_rate, lower, upper, ...) {
    standardGeneric("expected_cdo_loss")
  })

#' @rdname CalibrationParam-class
#' @aliases expected_cdo_loss,CalibrationParam-method
#'
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
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
  function(object, times, recovery_rate, lower, upper, ...,
      method = c("default", "CalibrationParam")) {
    method <- match.arg(method)
    qassert(recovery_rate, "N1[0,1]")
    qassert(lower, "N1[0,1]")
    assert_numeric(upper, lower = lower, upper = 1)
    g <- function(k, recovery_rate, lower, upper) {
      pmin(pmax((1 - recovery_rate) * k - lower, 0), upper - lower)
    }
    expected_value(object, times, g,
      recovery_rate = recovery_rate, lower = lower, upper = upper, ...)
  })

#' @describeIn CalibrationParam-class
#'   returns the expected portfolio CDS fair-value equation.
#'
#' @param discount_factors Discount factors for `times`
#' @param coupon Running coupon der CDO tranche.
#' @param upfront Upfront payment der CDO tranche.
#'
#' @export
setGeneric("expected_pcds_equation",
  function(object, times, discount_factors, recovery_rate, coupon, upfront, ...) {
    standardGeneric("expected_pcds_equation")
  })

#' @rdname CalibrationParam-class
#' @aliases expected_pcds_equation,CalibrationParam-method
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
function(object, times, discount_factors, recovery_rate, coupon, upfront, ...) {
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

#' @describeIn CalibrationParam-class
#'   returns the expected CDO fair-value equation.
#'
#' @export
setGeneric("expected_cdo_equation",
  function(object, times, discount_factors, recovery_rate, lower, upper, coupon, upfront, ...) {
    standardGeneric("expected_cdo_equation")
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
function(object, times, discount_factors, recovery_rate, coupon, upfront, ...) {
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
