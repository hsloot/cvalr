#' @include s4-ExtMO2FParam.R checkmate.R
NULL

#' Two-factor extendible Gaussian calibration parameters
#'
#' @description
#' [CalibrationParam-class]-class with two parameters for the extendible
#' equi-correlation Gaussian-copula model with Exponential  margins*(average)
#' default counting process*.
#'
#' @slot lambda A non-negative number for the marginal rate.
#' @slot nu A numeric number for the model specific dependence parameter
#'   (Pearson correlation between zero and one, use `rho` or `tau` to set
#'   dependence parameter).
#'
#' @details
#' The model is defined by the assumption that the multivariate default times
#' \eqn{\tau = (\tau_1, \ldots, \tau_d)} are extendible equi-correlation
#' Gaussian with Exponential margins.
#' The (internal) dependence parameter `nu` (*Pearson correlation*) has a
#' one-to-one relationship and can be replaced by *Spearman's Rho* `rho` or
#' *Kendall' Tau* `tau`. The possible range for `rho` or `tau` is from zero to
#' one (boundaries might not be included).
#' The link between the Pearson correlation coefficient and Spearman's Rho and
#' Kendall's Tau is
#' \itemize{
#'   \item \eqn{\rho = 2 \sin(\rho_S \cdot \pi / 6)} and
#'     \eqn{\rho_S = 6 / \pi \cdot \arcsin(\rho/2)}
#'   \item \eqn{\rho = \sin(\tau \cdot \pi / 2)} and
#'     \eqn{\tau = 2 / \pi \cdot \arcsin(\rho)}
#' }
#'
#' @export ExtGaussian2FParam
ExtGaussian2FParam <- setClass("ExtGaussian2FParam", # nolint
  contains = "CalibrationParam",
  slots = c(lambda = "numeric", nu = "numeric"))


setMethod("getLambda", "ExtGaussian2FParam",
  function(object) {
    object@lambda
  })
#' @importFrom checkmate qassert
setReplaceMethod("setLambda", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1(0,)")
    object@lambda <- value

    invisible(object)
  })

setMethod("getNu", "ExtGaussian2FParam",
  function(object) {
    object@nu
  })
#' @importFrom checkmate qassert
setReplaceMethod("setNu", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    object@nu <- value

    invisible(object)
  })

#' @importFrom checkmate qassert
setMethod("invRho", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    2 * sin(value * pi / 6)
  })

#' @importFrom checkmate qassert
setMethod("invTau", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    sin(value * pi / 2)
  })

setMethod("getRho", "ExtGaussian2FParam",
  function(object) {
    (6 / pi) * asin(getNu(object) / 2)
  })
#' @importFrom checkmate qassert
setReplaceMethod("setRho", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invRho(object, value)

    invisible(object)
  })

setMethod("getTau", "ExtGaussian2FParam",
  function(object) {
    (2 / pi) * asin(getNu(object))
  })
#' @importFrom checkmate qassert
setReplaceMethod("setTau", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    setNu(object) <- invTau(object, value)

    invisible(object)
  })


#' @importFrom checkmate qassert
setValidity("ExtGaussian2FParam",
  function(object) {
    qassert(object@lambda, "N1(0,)")
    qassert(object@nu, "N1[0,1]")

    invisible(TRUE)
  })


#' @describeIn ExtGaussian2FParam-class Constructor
#' @aliases initialize,ExtGaussian2FParam-method
#' @aliases initialize,ExtGaussian2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param dim Dimension.
#' @param lambda Marginal intensity.
#' @param nu (Internal) dependence parameter.
#' @param rho Bivariate Spearman's Rho.
#' @param tau Bivariate Kendall's Tau.
#'
#' @examples
#' ExtGaussian2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
setMethod("initialize", signature = "ExtGaussian2FParam",
  definition = function(.Object, # nolint
      dim, lambda, nu, rho, tau) {
    if (!missing(dim) && !missing(lambda) &&
          !(missing(nu) && missing(rho) && missing(tau))) {
      if (missing(nu)) {
        if (!missing(rho)) {
          nu <- invRho(.Object, rho)
        } else if (!missing(tau)) {
          nu <- invTau(.Object, tau)
        }
      }

      setDimension(.Object) <- dim
      setLambda(.Object) <- lambda
      setNu(.Object) <- nu
      validObject(.Object)
    }

    invisible(.Object)
  })



#' @describeIn ExtGaussian2FParam-class
#'    simulates the vector of *default times* and returns a matrix `x` with
#'    `dim(x) == c(n_sim, getDimension(object))`.
#' @aliases simulate_dt,ExtGaussian2FParam-method
#'
#' @inheritParams simulate_dt
#' @param method Simulation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#' @param n_sim Number of samples.
#'
#' @section Simulation:
#' The default times are sampled in a two-stage procedure: First a sample is
#' drawn from the equi-correlation Gaussian copula, see
#' [copula::normalCopula-class] and [copula::rCopula()]; then the results are
#' transformed using [stats::qexp()].
#'
#' @examples
#' parm <- ExtGaussian2FParam(5L, 8e-2, rho = 4e-1)
#' simulate_dt(parm, n_sim = 5L)
#'
#' @importFrom stats qexp
#' @importFrom copula normalCopula rCopula
#' @include utils.R
#'
#' @export
setMethod("simulate_dt", "ExtGaussian2FParam",
  function(object, ...,
      method = c("default", "ExtGaussian2FParam"), n_sim = 1e4L) {
    method <- match.arg(method)
    normalCopula(param = object@nu, dim = object@dim, dispstr = "ex") %>%
      rCopula(n_sim, .) %>%
      qexp(rate = object@lambda, lower.tail = FALSE)
  })

#' @describeIn ExtGaussian2FParam-class
#'   calculates the *probability vector* for the *average default count process*
#'   and returns a matrix `x` with `dim(x) == c(getDimension(object)+1L, length(times))`.
#' @aliases probability_distribution,ExtGaussian2FParam-method
#'
#' @inheritParams probability_distribution
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#'
#' @section Probability distribution:
#' The probability of \eqn{j} portfolio items being defaulted at time \eqn{t} is
#' calculated by numeric integration with [stats::integrate()] using
#' \deqn{
#'   \mathbb{P}\left( L_t = \frac{j}{d} \right)
#'     = \int_{-\infty}^{\infty} \binom{d}{j} G(t; x)^{j} [1 - G(t; x)]^{d-k}  \phi(x)\mathrm{d}x ,
#' }
#' with
#' \deqn{
#'   G(x; t)
#'     = \Phi\left( \frac{\Phi^{-1}(F(t)) - \sqrt{\nu}x}{\sqrt{1 - \nu}} \right)
#' }
#' and \eqn{F} being the Exponential distribution function for rate \eqn{\lambda}.
#'
#' @examples
#' parm <- ExtGaussian2FParam(5L, 8e-2, rho = 4e-1)
#' probability_distribution(parm, 0.3)
#' probability_distribution(parm, seq(0, 1, by = 0.25))
#' probability_distribution(parm, seq(0, 1, by = 0.25),
#'   method = "CalibrationParam", seed = 1623, sim_args(n_sim  = 1e2L))
#'
#' @importFrom stats integrate pexp pnorm dnorm qnorm
#' @importFrom checkmate qassert
#' @importFrom purrr map_dbl cross2
#' @include utils.R
#' @export
setMethod("probability_distribution", "ExtGaussian2FParam",
  function(object, times, ...,
      method = c("default", "ExtGaussian2FParam", "CalibrationParam")) {
    method <- match.arg(method)
    if (isTRUE("default" == method || "ExtGaussian2FParam" == method)) {
      qassert(times, "N+[0,)")
      d <- getDimension(object)
      lambda <- getLambda(object)
      nu <- getNu(object)
      out <- qnorm(pexp(times, rate = lambda)) %>%
        cross2(0L:d, .) %>%
        map_dbl(~{
          k <- .[[1]]
          z <- .[[2]]
          if (-Inf == z) {
            if (0L == k) {
              out <- 1
            } else {
              out <- 0
            }
          } else {
            lg <- function(x, lower_tail = TRUE) {
              pnorm((z - sqrt(nu) * x) / sqrt(1 - nu), log.p = TRUE, lower.tail = lower_tail)
            }
            int_fn <- function(x) {
              v_multiply_binomial_coefficient(
                exp(k * lg(x, lower_tail = TRUE) + (d - k) * lg(x, lower_tail = FALSE)), d, k) *
                dnorm(x)
            }
            int_res <- integrate(
              int_fn, lower = -Inf, upper = Inf, rel.tol = .Machine$double.eps^0.5)
            out <- int_res$value
          }
          out
        }) %>%
        matrix(nrow = d+1L, ncol = length(times))
    } else {
      out <- callNextMethod(object, times, ..., method = method)
    }

    out
  })

#' @describeIn ExtGaussian2FParam-class
#'   returns the expected portfolio CDS loss for a specific time-point.
#' @aliases expected_pcds_loss,ExtGaussian2FParam-method
#'
#' @inheritParams probability_distribution
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#'
#' @inheritSection ExtMO2FParam-class Expected portfolio CDS loss
#'
#' @examples
#' parm <- ExtGaussian2FParam(75L, 8e-2, rho = 4e-1)
#' expected_pcds_loss(parm, times = 0.25, recovery_rate = 0.4)
#' expected_pcds_loss(parm, times = seq(0, 1, by = 0.25), recovery_rate = 0.4)
#' expected_pcds_loss(parm, times = seq(0, 1, by = 0.25), recovery_rate = 0.4,
#'   method = "CalibrationParam")
#' expected_pcds_loss(parm, times = seq(0, 1, by = 0.25), recovery_rate = 0.4,
#'   method = "CalibrationParam",
#'   pd_args = list(method = "CalibrationParam", seed = 1623,
#'     sim_args = list(n_sim = 1e2L)))
#'
#' @importFrom stats pexp
#' @importFrom checkmate qassert
#'
#' @export
setMethod("expected_pcds_loss", "ExtGaussian2FParam",
  function(object, times, recovery_rate, ...,
      method = c("default", "ExtGaussian2FParam", "CalibrationParam")) {
    method <- match.arg(method)
    if (isTRUE("default" == method || "ExtGaussian2FParam" == method)) {
      qassert(times, "N+[0,)")
      qassert(recovery_rate, "N1[0,1]")
      out <- (1 - recovery_rate) * pexp(times, rate = object@lambda)
    } else {
      out <- callNextMethod(object, times, recovery_rate, ..., method = method)
    }

    out
  })

#' @describeIn ExtGaussian2FParam-class
#'   returns the expected CDO loss for a specific time-point.
#' @aliases expected_cdo_loss,ExtGaussian2FParam-method
#'
#' @inheritParams probability_distribution
#' @param method Calculation method (either `"default"` or the name of the
#'   class whose implementation should be used).
#'
#' @section Expected CDO tranche loss:
#' The *expected loss of a CDO tranche* with *upper and lower attachment points*
#' \eqn{l} and \eqn{u} and *recovery rate* \eqn{R} is calculated using that
#' \deqn{
#'   \mathbb{E}[g(L_t)]
#'     = (1 - R) \left( C_2\left(\left[1 - \frac{l}{1-R} \right] \vee 0, F(t); -\sqrt{1 - \nu}\right) - C_2\left(\left[1 - \frac{u}{1-R} \right] \vee 0, F(t); -\sqrt{1 - \nu}\right) \right)
#' }
#' with \eqn{g(x) = \{[(1 - R) x - l] \vee (u-l)\} \wedge 1}, \eqn{C_2} being
#' the bivariate Gaussian copula, and \eqn{F} being the Exponential distribution
#' function for \eqn{\lambda}.
#'
#' @examples
#' parm <- ExtGaussian2FParam(5L, 8e-2, rho = 4e-1)
#' expected_cdo_loss(parm, times = 0.25, recovery_rate = 0.4, lower = 0.1, upper = 0.2)
#' expected_cdo_loss(parm, times = seq(0, 1, by = 0.25),
#'   recovery_rate = 0.4, lower = 0.1, upper = 0.2)
#' expected_cdo_loss(parm, times = seq(0, 1, by = 0.25),
#'   recovery_rate = 0.4, lower = 0.1, upper = 0.2, method = "CalibrationParam")
#' expected_cdo_loss(parm, times = seq(0, 1, by = 0.25),
#'   recovery_rate = 0.4, lower = 0.1, upper = 0.2, method = "CalibrationParam",
#'   pd_args = list(method = "CalibrationParam", seed = 1623,
#'     sim_args = list(n_sim = 1e2L)))
#'
#' @importFrom copula normalCopula pCopula
#' @importFrom stats pexp
#' @importFrom checkmate qassert assert_numeric
#'
#' @export
setMethod("expected_cdo_loss", "ExtGaussian2FParam",
  function(
      object, times, recovery_rate, lower, upper, ...,
      method = c("default", "ExtGaussian2FParam", "CalibrationParam")) {
    method <- match.arg(method)
    if (isTRUE("default" == method || "ExtGaussian2FParam" == method)) {
      qassert(times, "N+[0,)")
      qassert(recovery_rate, "N1[0,1]")
      qassert(lower, "N1[0,1]")
      assert_numeric(upper, lower = lower, upper = 1)

      cop <- normalCopula(-sqrt(1-object@nu))
      out <- pexp(times, rate = object@lambda) %>%
        map_dbl(~{
          u_left <- c(pmax(1 - lower / (1 - recovery_rate), 0), .)
          u_right <- c(pmax(1 - upper / (1 - recovery_rate), 0), .)
          (1 - recovery_rate) * (pCopula(u_left, cop) - pCopula(u_right, cop))
        })
    } else {
      out <- callNextMethod(object, times, recovery_rate, lower, upper, ...)
    }

    out
  })


#' @describeIn ExtMOParam-class Display the object.
#' @aliases show,ExtGaussian2FParam-method
#'
#' @inheritParams methods::show
#'
#' @export
setMethod("show", "ExtGaussian2FParam",
 function(object) {
   cat(sprintf("An object of class %s\n", classLabel(class(object))))
   cat(sprintf("Dimension: %i\n", getDimension(object)))
   cat(sprintf(
     "Lambda: %s, Rho: %s, Tau: %s\n",
     format(getLambda(object)), format(getRho(object)), format(getTau(object))))
  })
