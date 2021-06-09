#' @importFrom methods setClass setValidity setGeneric setMethod setReplaceMethod
#'   validObject is as new callNextMethod show classLabel hasMethod
#' @include checkmate.R
NULL

# nolint start
ERR_MSG_DIM <- "`dim` must be scalar integer >= 2"
ERR_MSG_REQ_IMPL <- "`%s` requires implementation of `%s`"
MSG_IGN_PARM <- "parameter `%s` is ignored for `%s`"
WRN_MISSING_PARM <- "parameter `%s` missing; using default %s"
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
#'   [JoeExtArch2FParam-class]
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
#' @param times A non-negative numeric vector of timepoints.
#'
#' @export
setGeneric("simulate_adcp",
  function(object, times, ...) {
    standardGeneric("simulate_adcp")
  })

#' @describeIn CalibrationParam-class
#'   simulates the *average default counting process* and returns a
#'   matrix `x` with `dim(x) == c(n_sim, length(times))`.
#' @aliases simulate_adcp,CalibrationParam-methods
#'
#' @examples
#' parm <- ExMarkovParam(rmo::exQMatrix(rmo::AlphaStableBernsteinFunction(0.4), d = 5L))
#' simulate_adcp(parm, 1, n_sim = 10L)
#' simulate_adcp(parm, seq(0, 5, by = 0.25), n_sim = 10L)
#'
#' @importFrom stats rexp
#' @include RcppExports.R
#'
#' @export
setMethod("simulate_adcp", "CalibrationParam",
  function(object, times, ...) {
    qassert(times, "N+[0,)")

    Rcpp__dt2adcp(simulate_dt(object, ...), times)
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

#' @describeIn CalibrationParam-class
#'   calculates the *expected value* for the loss based on the *average default
#'   count process* for given timepoints and returns a vector `x` with
#'   `length(x) == length(times)`.
#'
#' @section Expected value:
#' The *expectated value* of finite linear transformations of the *average default counting process*
#' \eqn{L} under a transformation and a linear aggregation can be calculated with
#' [expected_value()]; i.e. the value
#' \deqn{
#'   \mathbb{E}[diag(A^\top \cdot g(L_t))] , \quad t = (t_1, \ldots, t_m) \geq 0.
#' }
#' with \eqn{g : \mathbb{R} \to \mathbb{R}^p} \eqn{A = (a_1, \ldots, a_p) \in R^{m \times p}}.
#'
#' * For classes with an implementation of `probability_distribution`, `method == "prob"` can be
#'   used. In this case the transformation `g` (generalized such that \eqn{g(\mathbb{R}^m) \in
#'   \mathbb{R}^{m \times p}}) should be provided via `.trans_v(x)` which takes an
#'   \eqn{\mathbb{R}^m} vector and returns a \eqn{\mathbb{R}^{m\times p}} matrix. Furthermore, the
#'   (diagonal of the) linear aggregation should be provided via `.lagg_ev(x)` which takes a
#'   \eqn{\mathbb{R}^{m \times p}} matrix and returns an \eqn{\mathbb{R}^p} vector. Both functions
#'   should be vectorized appropriately.
#' * For all classes, `method == "prob"` can be used. In this case a function `.simulate_pv` should
#'   be provided. This function should draw samples of \eqn{L_t^{(i)}} and  return the values of
#'   \eqn{(a_j^\top \cdot g_j(L_t^{(i)}))} in a \eqn{\mathbb{R}^{n \times p}} matrix. A named list
#'   of functions can be provided via `attrs` those should work simular to [mean()]; their results
#'   will be bound to the return values as attributes.
#'
#' The methods [expected_pcds_equation()] and [expected_cdo_equation()] are convenience-wrappers for
#' expected *portfolio CDS* and *CDO* payoff-equations.
#'
#' @export
setGeneric("expected_value",
  function(object, times, ...) {
    standardGeneric("expected_value")
  })

#' @describeIn CalibrationParam-class
#'   calculates the *expected value* for the loss based on the *average default
#'   count process* for given timepoints and returns a vector `x` with
#'   `length(x) == length(times)`.
#' @aliases expected_value,CalibrationParam-method
#'
#' @param method Calculation method (either `"default"`, `"prob"` (requires implementation of
#'   `probability_distribution`), or `"mc"`).
#' @param n_sim Number of samples.
#' @param attrs A named list with functions which are applied to a vector of *present values*.
#' @param .trans_v Internal parameter, not independent for the user.
#' @param .lagg_ev Internal parameter, not independent for the user.
#' @param .simulate_pv Internal parameter, not independent for the user.
#'
#' @examples
#' expected_value(CuadrasAugeExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.25,
#'   .trans_v = function(x) x, .lagg_ev = function(x) x)
#' expected_value(AlphaStableExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), seq(0, 5, by = 0.25),
#'   .trans_v = function(x) x, .lagg_ev = function(x) apply(x, 2L, mean))
#' expected_value(PoissonExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.25,
#'   .trans_v = function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2),
#'   .lagg_ev = function(x) apply(x, 2L, mean))
#' expected_value(ExponentialExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), seq(0, 5, by = 0.25),
#'   .trans_v = function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2),
#'   .lagg_ev = function(x) apply(x, 2L, mean))
#'
#' expected_value(CuadrasAugeExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), 0.25,
#'   method = "mc", n_sim = 1e4L,
#'   .simulate_pv = function(object, times, n_sim) simulate_adcp(object, times, n_sim = n_sim))
#' expected_value(CuadrasAugeExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), seq(0, 5, by = 0.25),
#'   method = "mc", n_sim = 1e4L,
#'   .simulate_pv = function(object, times, n_sim) simulate_adcp(object, times, n_sim = n_sim))
#' expected_value(CuadrasAugeExtMO2FParam(
#'   dim = 50L, lambda = 0.05, rho = 0.4), seq(0, 5, by = 0.25),
#'   function(x) pmin(pmax(0.6 * x - 0.1, 0), 0.2),
#'   method = "mc", n_sim = 1e4L,
#'   .simulate_pv = function(object, times, n_sim) simulate_adcp(object, times, n_sim = n_sim),
#'   attrs = list(sd = function(x) sd(x) / sqrt(length(x))))
#'
#' @importFrom checkmate qassert assert_function assert_list
#' @importFrom purrr map_dbl
#' @include utils-math.R
#' @export
setMethod("expected_value", "CalibrationParam", # nolint
  function(object, times, ...,
      method = c("default", "prob", "mc"), n_sim = 1e4L, attrs = NULL,
      .trans_v = NULL, .lagg_ev = NULL, .simulate_pv = NULL) {
    method <- match.arg(method)
    qassert(times, "N+[0,)")

    if ("default" == method) {
      if (hasMethod("probability_distribution", class(object))) {
        method <- "prob"
      } else {
        method <- "mc"
      }
    }

    if ("prob" == method) {
      assert_function(.trans_v, null.ok = FALSE, args = "x", ordered = TRUE)
      assert_function(.lagg_ev, null.ok = FALSE, args = "x", ordered = TRUE)
      assert_list(attrs, types = "function", names = "named", null.ok = TRUE)
      if (!hasMethod("probability_distribution", class(object))) {
        stop(sprintf(
          ERR_MSG_REQ_IMPL, "expected_value(., method == \"prob\")", "probability_distribution"))
      }
      if (!missing(n_sim) && isTRUE(options("cvalr.enable_messages"))) {
        stop(sprintf(MSG_IGN_PARM, "n_sim", "expected_value(., method = \"prob\")"))
      }
      if (!missing(attrs) && isTRUE(options("cvalr.enable_messages"))) {
        stop(sprintf(MSG_IGN_PARM, "attrs", "expected_value(., method = \"prob\")"))
      }

      d <- getDimension(object)
      vals <- .trans_v((0:d) / d)
      probs <- probability_distribution(object, times)

      out <- .lagg_ev(dot(probs, vals, simplify = FALSE))
    } else {
      qassert(n_sim, "X1[1,)")
      assert_function(.simulate_pv, null.ok = FALSE, args = c("object", "times", "n_sim"),
        ordered = TRUE)
      if (missing(n_sim) && isTRUE(options("cvalr.enable_warnings"))) {
        warning(sprintf(WRN_MISSING_PARM, "n_sim", format(n_sim)))
      }
      pv <- .simulate_pv(object, times, n_sim = n_sim)
      out <- apply(pv, 2L, mean)
      if (!is.null(attrs)) {
        for (i in seq_along(attrs)) {
          label <- names(attrs)[[i]]
          fn <- attrs[[i]]
          attr(out, label) <- apply(pv, 2L, fn)
        }
      }
    }

    out
  })

#' @describeIn CalibrationParam-class
#'   calculates the *payoff equation* for a *portfolio CDS* (vectorized w.r.t.
#'   the argumentes `recovery_rate`, `coupon`, and `upfront`).
#'
#' @param discount_factors Non-negative numeric vector for the discount factors for the timepoints.
#' @param recovery_rate Non-negative number between zero and one for the *recovery rate*..
#' @param coupon Numeric number for the running coupon.
#' @param upfront Numeric number for the upfront payment.
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
#'   AlphaStableExtMO2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   times = seq(0, 5, by = 0.25), discount_factors = rep(1, 21),
#'   recovery_rate = 0.4, coupon = 0.08, upfront = 0)
#' expected_pcds_equation(
#'   AlphaStableExtMO2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   times = seq(0, 5, by = 0.25), discount_factors = rep(1, 21),
#'   recovery_rate = 0.4, coupon = rep(0.08, 4), upfront = rep(0, 4))
#' expected_pcds_equation(
#'   AlphaStableExtMO2FParam(dim = 75, lambda = 0.05, rho = 0.6),
#'   times = seq(0, 5, by = 0.25), discount_factors = rep(1, 21),
#'   recovery_rate = rep(0.4, 4), coupon = rep(0.08, 4), upfront = rep(0, 4))
#'
#' @importFrom checkmate qassert
#' @importFrom purrr pmap_dbl
#' @importFrom vctrs vec_size_common vec_recycle
#'
#' @export
setMethod("expected_pcds_equation", "CalibrationParam",
function(object, times, discount_factors, recovery_rate, coupon, upfront, ...) {
  qassert(times, "N+[0,)")
  qassert(discount_factors, paste0("N", length(times), "[0,)"))
  qassert(recovery_rate, "N+[0,1]")
  qassert(coupon, "N+")
  qassert(upfront, "N+")

  p <- vec_size_common(recovery_rate, coupon, upfront)
  recovery_rate <- vec_recycle(recovery_rate, p)
  coupon <- vec_recycle(coupon, p)
  upfront <- vec_recycle(upfront, p)

  .trans_v <- function(x) {
    Rcpp__trans_v_pcds(x, recovery_rate)
  }
  .lagg_ev <- function(x) {
    Rcpp__lagg_ev_pcds(x, times, discount_factors, recovery_rate, coupon, upfront)
  }
  .simulate_pv <- function(object, times, n_sim) {
    simulate_adcp(object, times, n_sim = n_sim) %>%
      Rcpp__adcp2peqpv_pcds(times, discount_factors, recovery_rate, coupon, upfront)
  }

  expected_value(
    object, times, ..., .trans_v = .trans_v, .lagg_ev = .lagg_ev, .simulate_pv = .simulate_pv)
  })

#' @describeIn CalibrationParam-class
#'   calculates the *payoff equation* for a *CDO* (vectorized w.r.t. the
#'   argumentes `recovery_rate`, `coupon`, and `upfront`).
#'
#' @param lower Non-negative number between zero and one for the *lower attachment point*.
#' @param upper Non-negative number between zero and one for the *upper attachment point*.
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
#'   times = seq(0, 5, by = 0.25), discount_factors = rep(1, 21),
#'   recovery_rate = 0.4, lower = c(0, 0.1, 0.2, 0.35),
#'   upper = c(0.1, 0.2, 0.35, 1), coupon = 0.08, upfront = 0
#' )
#'
#' @importFrom checkmate qassert
#' @importFrom purrr pmap_dbl
#' @importFrom vctrs vec_size_common vec_recycle
#'
#' @export
setMethod("expected_cdo_equation", "CalibrationParam",
function(object, times, discount_factors, recovery_rate, lower, upper, coupon, upfront, ...) {
  qassert(times, "N+[0,)")
  qassert(discount_factors, paste0("N", length(times), "[0,)"))
  qassert(recovery_rate, "N+[0,1]")
  qassert(lower, "N+[0,1]")
  qassert(upper, "N+[0,1]")
  qassert(upper - lower, "N+[0,1]")
  qassert(coupon, "N+")
  qassert(upfront, "N+")

  p <- vec_size_common(recovery_rate, lower, upper, coupon, upfront)
  recovery_rate <- vec_recycle(recovery_rate, p)
  lower <- vec_recycle(lower, p)
  upper <- vec_recycle(upper, p)
  coupon <- vec_recycle(coupon, p)
  upfront <- vec_recycle(upfront, p)

  .trans_v <- function(x) {
    Rcpp__trans_v_cdo(x, recovery_rate, lower, upper)
  }
  .lagg_ev <- function(x) {
    Rcpp__lagg_ev_cdo(x, times, discount_factors, recovery_rate, lower, upper, coupon, upfront)
  }
  .simulate_pv <- function(object, times, n_sim) {
    simulate_adcp(object, times, n_sim = n_sim) %>%
      Rcpp__adcp2peqpv_cdo(times, discount_factors, recovery_rate, lower, upper, coupon, upfront)
  }

  expected_value(
    object, times, ..., .trans_v = .trans_v, .lagg_ev = .lagg_ev, .simulate_pv = .simulate_pv)
  })
