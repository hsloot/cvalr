#' @importFrom methods setClass setValidity setGeneric setMethod setReplaceMethod
#'   validObject is as new callNextMethod show classLabel
#' @include checkmate.R
NULL

# nolint start
ERR_MSG_DIM <- "`dim` must be scalar integer >= 2"
# nolint end

#' Virtual super-class for calibration parameters
#'
#' @description
#' `CalibrationParam` provides a simple interface to calculate *expected values*
#' and *pricing equations* for *portfolio CDS'* and *CDO's* with different
#' models.
#'
#' @slot dim The dimension (number of portfolio items).
#'
#' @seealso [ExMarkovParam-class] [ExMOParam-class] [ExtMOParam-class]
#'   [CuadrasAugeExtMO2FParam-class] [AlphaStableExtMO2FParam-class]
#'   [PoissonExtMO2FParam-class] [ExponentialExtMO2FParam-class]
#'   [ExtArch2FParam-class] [ClaytonExtArch2FParam-class]
#'   [FrankExtArch2FParam-class] [GumbelExtArch2FParam-class]
#'   [JoeExtArch2FParam-class] [expected_pcds_loss()], [expected_cdo_loss()]
#'
#' @export
setClass("CalibrationParam", # nolint
  contains = "VIRTUAL",
  slots = c(dim = "integer"))


setGeneric("getDimension",
  function(object) {
    standardGeneric("getDimension")
  })

setMethod("getDimension", "CalibrationParam",
  function(object) {
    object@dim
  })

setGeneric("setDimension<-",
  function(object, value) {
    standardGeneric("setDimension<-")
  })

#' @importFrom checkmate qassert
setReplaceMethod("setDimension", "CalibrationParam",
  function(object, value) {
    qassert(value, "X1[2,)")
    object@dim <- as.integer(value)

    invisible(object)
  })


#' @importFrom checkmate qtest
setValidity("CalibrationParam",
  function(object) {
    if (!qtest(object@dim, "I1[2,)")) {
      return(ERR_MSG_DIM)
    }

    invisible(TRUE)
  })


#' @describeIn CalibrationParam-class
#'    simulates the vector of *default times* and returns a matrix `x` with
#'    `dim(x) == c(n_sim, getDimension(object))`.
#'
#' @param object A [CalibrationParam-class]-object.
#' @param ... Pass-through parameters.
#'
#' @export
setGeneric("simulate_dt",
  function(object, ...) {
    standardGeneric("simulate_dt")
  })

#' @describeIn CalibrationParam-class
#'   simulates the *average default counting process* and returns a
#'   matrix `x` with `dim(x) == c(n_sim, length(times))`.
#'
#' @param object A [CalibrationParam-class]-object.
#' @param times A non-negative numeric vector of timepoints.
#' @param ... Pass-through parameters.
#'
#' @export
setGeneric("simulate_adcp",
  function(object, times, ...) {
    standardGeneric("simulate_adcp")
  })

#' @describeIn CalibrationParam-class
#'   simulates the *average default counting process* and returns a
#'   matrix `x` with `dim(x) == c(n_sim, length(times))`.
#' @aliases simulate_adcp,CalibrationParam-method
#'
#' @examples
#' parm <- ExMarkovParam(rmo::exQMatrix(rmo::AlphaStableBernsteinFunction(0.4), d = 5L))
#' simulate_adcp(parm, 1, n_sim = 10L)
#' simulate_adcp(parm, seq(0, 5, by = 0.25), n_sim = 10L)
#'
#' @importFrom stats rexp
#' @include utils.R
#'
#' @export
setMethod("simulate_adcp", "CalibrationParam",
  function(object, times, ...) {
     simulate_dt(object, ...) %>%
      dt2adcp(times)
  })

#' @describeIn CalibrationParam-class
#'   calculates the *probability vector* for the *average default count process*
#'   and returns a matrix `x` with `dim(x) == c(getDimension(object)+1L, length(times))`.
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
#'   calculates the *probability vector* for the *average default count process*
#'   and returns a matrix `x` with `dim(x) == c(getDimension(object)+1L, length(times))`.
#' @aliases probability_distribution,CalibrationParam-method
#'
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param seed Numeric number (if not NULL, it is used to set the seed prior to
#'   Monte-Carlo estimation of the probability distribution).
#' @param sim_args List with pass-through parameters for [simulate_adcp()].
#'
#' @importFrom checkmate qassert assert_number assert_list
#' @importFrom rlang exec
#' @include utils.R
#'
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
    exec(simulate_adcp, object = object, times = times, !!!sim_args) %>%
      adcp2epd(getDimension(object))
  })

#' @describeIn CalibrationParam-class
#'   calculates the *expected value* for the loss based on the *average default
#'   count process* for given timepoints and returns a vector `x` with
#'   `length(x) == length(times)`.
#'
#' @param object A [CalibrationParam-class]-object.
#' @param g A payoff transformation function.
#' @param ... Pass-through parameter.
#'
#' @section Expected value:
#' The *expectated value* of the *average default counting process* \eqn{L}
#' under a payoff-transformation can be calculated with [expected_value()]; i.e. the
#' value
#' \deqn{
#'   \mathbb{E}[g(L_t)] , \quad t \geq 0.
#' }
#' For a *portfolio CDS* choose \eqn{g(x) = (1 - R) x} and for a *CDO* tranche with
#' attachment points \eqn{l < u} and choose
#' \eqn{g(x) = min{\{ \max{\{ (1 - R) x - l, 0 \}}, u - l \}}}, where \eqn{R} is
#' the recovery rate.
#'
#' * [expected_value()] can be used for general payoff-transformations of the
#'   *average default counting process*.
#' * [expected_pcds_loss()] and [expected_cdo_loss()] are convenience-wrappers
#'   for *portfolio CDS* and *CDO* payoff-transformations.
#'
#' @export
setGeneric("expected_value",
  function(object, times, g, ...) {
    standardGeneric("expected_value")
  })

#' @describeIn CalibrationParam-class
#'   calculates the *expected value* for the loss based on the *average default
#'   count process* for given timepoints and returns a vector `x` with
#'   `length(x) == length(times)`.
#' @aliases expected_value,CalibrationParam-method
#'
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param pd_args A list of parameters for [probability_distribution()].
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
#' @importFrom purrr map_dbl
#' @importFrom rlang exec
#' @include utils-math.R
#' @export
setMethod("expected_value", "CalibrationParam",
  function(object, times, g, ...,
      method = c("default", "CalibrationParam"), pd_args = NULL) {
    method <- match.arg(method)
    qassert(times, "N+[0,)")
    assert_function(g)
    assert_list(pd_args, null.ok = TRUE)
    d <- getDimension(object)
    vals <- map_dbl((0:d) / d, g, ...)
    probs <- exec(probability_distribution, object = object, times = times, !!!pd_args)

    dot(vals, probs)
  })

#' @describeIn CalibrationParam-class
#'   calculates the *expected value* for the *portfolio CDS loss* based on the
#'   *average default count process* for given timepoints and returns a vector
#'   `x` with `length(x) == length(times)`.
#'
#' @param recovery_rate Non-negative number between zero and one for the recovery rate
#'   of the portfolio CDS/CDO.
#'
#' @export
setGeneric("expected_pcds_loss",
  function(object, times, recovery_rate, ...) {
    standardGeneric("expected_pcds_loss")
  })

#' @describeIn CalibrationParam-class
#'   calculates the *expected value* for the *portfolio CDS loss* based on the
#'   *average default count process* for given timepoints and returns a vector
#'   `x` with `length(x) == length(times)`.
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
#'
#' @export
setMethod("expected_pcds_loss", "CalibrationParam",
  function(object, times, recovery_rate, ...,
      method = c("default", "CalibrationParam")) {
    method <- match.arg(method)
    qassert(recovery_rate, "N1[0,1]")
    g_pcds <- function(x, recovery_rate) {
      (1 - recovery_rate) * x
    }

    expected_value(object, times, g_pcds, recovery_rate = recovery_rate, ...)
  })

#' @describeIn CalibrationParam-class
#'   calculates the *expected value* for the *CDO loss* based on the *average
#'   default count process* for given timepoints and returns a vector `x` with
#'   `length(x) == length(times)`.
#'
#' @param lower Non-negative number between zero and one for the lower attachment point
#'   of the CDO tranche.
#' @param upper Non-negative number between `lower` and one for the upper attachment point
#'   of the CDO tranche.
#'
#' @export
setGeneric("expected_cdo_loss",
  function(object, times, recovery_rate, lower, upper, ...) {
    standardGeneric("expected_cdo_loss")
  })

#' @describeIn CalibrationParam-class
#'   calculates the *expected value* for the *CDO loss* based on the *average
#'   default count process* for given timepoints and returns a vector `x` with
#'   `length(x) == length(times)`.
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
#'
#' @export
setMethod("expected_cdo_loss", "CalibrationParam",
  function(object, times, recovery_rate, lower, upper, ...,
      method = c("default", "CalibrationParam")) {
    method <- match.arg(method)
    qassert(recovery_rate, "N1[0,1]")
    qassert(lower, "N1[0,1]")
    assert_numeric(upper, lower = lower, upper = 1)
    g_cdo <- function(x, recovery_rate, lower, upper) {
      pmin(pmax((1 - recovery_rate) * x - lower, 0), upper - lower)
    }
    expected_value(
      object, times, g_cdo, recovery_rate = recovery_rate, lower = lower, upper = upper, ...)
  })

#' @describeIn CalibrationParam-class
#'   calculates the *payoff equation* for a *portfolio CDS* (vectorized w.r.t.
#'   the argumentes `recovery_rate`, `coupon`, and `upfront`).
#'
#' @param discount_factors Non-negative numeric vector for the discount factors for the timepoints.
#' @param coupon Numeric number for the running coupon of the portfolio CDS / CDO tranche.
#' @param upfront Numeric number for the upfront payment of the portfolio CDS / CDO tranche.
#'
#' @export
setGeneric("expected_pcds_equation",
  function(object, times, discount_factors, recovery_rate, coupon, upfront, ...) {
    standardGeneric("expected_pcds_equation")
  })

#' @describeIn CalibrationParam-class
#'   calculates the *payoff equation* for a *portfolio CDS* (vectorized w.r.t.
#'   the argumentes `recovery_rate`, `coupon`, and `upfront`).
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
#' @importFrom checkmate qassert
#' @importFrom purrr pmap_dbl
#'
#' @export
setMethod("expected_pcds_equation", "CalibrationParam",
function(object, times, discount_factors, recovery_rate, coupon, upfront, ...) {
  qassert(coupon, "N+")
  qassert(upfront, "N+")
  pricing_equation <- function(
      parm, times, discount_factors, recovery_rate, coupon, upfront, ...) {
    portfolio_cds_equation(
      expected_pcds_loss(parm, times, recovery_rate, ...),
      times, discount_factors, recovery_rate, coupon, upfront)
  }

  pmap_dbl(
    list(recovery_rate = recovery_rate, coupon = coupon, upfront = upfront),
    pricing_equation,
    parm = object, times = times, discount_factors = discount_factors, ...)
  })

#' @describeIn CalibrationParam-class
#'   calculates the *payoff equation* for a *CDO* (vectorized w.r.t. the
#'   argumentes `recovery_rate`, `coupon`, and `upfront`).
#'
#' @export
setGeneric("expected_cdo_equation",
  function(object, times, discount_factors, recovery_rate, lower, upper, coupon, upfront, ...) {
    standardGeneric("expected_cdo_equation")
  })

#' @describeIn CalibrationParam-class
#'   calculates the *payoff equation* for a *CDO* (vectorized w.r.t. the
#'   argumentes `recovery_rate`, `coupon`, and `upfront`).
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
#' @importFrom checkmate qassert
#' @importFrom purrr pmap_dbl
#'
#' @export
setMethod("expected_cdo_equation", "CalibrationParam",
function(object, times, discount_factors, recovery_rate, lower, upper, coupon, upfront, ...) {
  qassert(coupon, "N+")
  qassert(upfront, "N+")
  pricing_equation <- function(
      parm, times, discount_factors, recovery_rate, lower, upper, coupon, upfront, ...) {
    cdo_equation(
      expected_cdo_loss(
        parm, times, recovery_rate, lower, upper, ...),
      times, discount_factors, lower, upper, coupon, upfront)
  }
  pmap_dbl(
    list(recovery_rate = recovery_rate, lower = lower, upper = upper,
         coupon = coupon, upfront = upfront),
    pricing_equation,
    parm = object, times = times, discount_factors = discount_factors, ...)
  })
