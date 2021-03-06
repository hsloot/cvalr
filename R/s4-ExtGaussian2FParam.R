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
    if (!qtest(object@lambda, "N1(0,)")) {
      return(ERR_MSG_LAMBDA)
    }
    if (!qtest(object@nu, "N1")) {
      return(sprintf(ERR_MSG_NU1_INTERVAL, "[0,1]"))
    }

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
#' ExtGaussian2FParam()
#' ExtGaussian2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
setMethod("initialize", signature = "ExtGaussian2FParam",
  definition = function(.Object, # nolint
      dim, lambda, nu, rho, tau) {
    if (!missing(dim) && !missing(lambda) &&
          !(missing(nu) && missing(rho) && missing(tau))) {
      setDimension(.Object) <- dim
      setLambda(.Object) <- lambda
      if (!missing(nu)) {
        setNu(.Object) <- nu
      } else if (!missing(rho)) {
        setRho(.Object) <- rho
      } else if (!missing(tau)) {
        setTau(.Object) <- tau
      }

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
#' @param n_sim Number of samples.
#'
#' @section Simulation:
#' The default times are sampled in a two-stage procedure: First a sample is drawn from the
#' equi-correlation Gaussian copula, see [copula::normalCopula-class] and [copula::rCopula()]; then
#' the results are transformed using [stats::qexp()].
#'
#' @examples
#' parm <- ExtGaussian2FParam(5L, 8e-2, rho = 4e-1)
#' simulate_dt(parm, n_sim = 5L)
#'
#' @importFrom stats qexp pnorm
#' @importFrom mvtnorm rmvnorm
#' @include utils.R
#'
#' @export
setMethod("simulate_dt", "ExtGaussian2FParam",
  function(object, ..., n_sim = 1e4L) {
    d <- getDimension(object)
    corr <- matrix(getNu(object), nrow = d, ncol = d)
    diag(corr) <- rep(1, d)
    qexp(pnorm(rmvnorm(n_sim, sigma = corr, checkSymmetry = FALSE)),
      rate = getLambda(object), lower.tail = FALSE)
  })

#' @describeIn ExtGaussian2FParam-class
#'   calculates the *probability vector* for the *average default count process*
#'   and returns a matrix `x` with `dim(x) == c(getDimension(object)+1L, length(times))`.
#' @aliases probability_distribution,ExtGaussian2FParam-method
#'
#' @inheritParams probability_distribution
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
  function(object, times, ...) {
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

    out
  })

#' @describeIn ExtGaussian2FParam-class
#'   calculates the *payoff equation* for a *portfolio CDS* (vectorized w.r.t.
#'   the argumentes `recovery_rate`, `coupon`, and `upfront`).
#' @aliases expected_pcds_equation,ExtGaussian2FParam-method
#'
#' @inheritParams expected_pcds_equation
#'
#' @inheritSection ExtMO2FParam-class Expected portfolio CDS loss
#'
#' @examples
#' parm <- ExtGaussian2FParam(75L, 8e-2, rho = 4e-1)
#' expected_pcds_equation(
#'   parm, times = seq(25e-2, 5, by = 25e-2), discount_factors = rep(1, 20L), recovery_rate = 0.4,
#'   coupon = 1e-1, upfront = 0)
#' expected_pcds_equation(
#'   parm, times = seq(25e-2, 5, by = 25e-2), discount_factors = rep(1, 20L), recovery_rate = 0.4,
#'   coupon = 1e-1, upfront = 0, method = "mc", n_sim = 1e1)
#'
#' @importFrom stats pexp
#' @importFrom checkmate assert check_numeric
#' @importFrom vctrs vec_size_common vec_recycle
#' @include RcppExports.R
#'
#' @export
setMethod("expected_pcds_equation", "ExtGaussian2FParam",
  function(object, times, discount_factors, recovery_rate, coupon, upfront, ...,
      method = c("default", "prob", "mc")) {
    method <- match.arg(method)
    if ("default" == method) {
      qassert(times, "N+[0,)")
      qassert(discount_factors, paste0("N", length(times), "[0,)"))
      qassert(recovery_rate, "N+[0,1]")
      qassert(coupon, "N+")
      qassert(upfront, "N+")

      p <- vec_size_common(recovery_rate, coupon, upfront)
      recovery_rate <- vec_recycle(recovery_rate, p)
      coupon <- vec_recycle(coupon, p)
      upfront <- vec_recycle(upfront, p)

      x <- pexp(times, rate = getLambda(object)) %*% t(1 - recovery_rate)

      out <- Rcpp__lagg_ev_pcds(x, times, discount_factors, recovery_rate, coupon, upfront)
    } else {
      out <- callNextMethod(
        object, times, discount_factors, recovery_rate, coupon, upfront, ..., method = method)
    }

    out
  })

#' @describeIn ExtGaussian2FParam-class
#'   calculates the *payoff equation* for a *CDO* (vectorized w.r.t. the
#'   argumentes `recovery_rate`, `coupon`, and `upfront`).
#' @aliases expected_cdo_equation,ExtGaussian2FParam-method
#'
#' @inheritParams expected_cdo_equation
#' @param method Calculation method (either `"default"`, `"prob"`, or `"mc"`).
#'
#' @section Expected CDO tranche loss:
#' The *expected loss of a CDO tranche* with *upper and lower attachment points*
#' \eqn{l} and \eqn{u} and *recovery rate* \eqn{R} is calculated using that
#' \deqn{
#'   \mathbb{E}[g(L_t)]
#'     = (1 - R) \left(
#'       C_2\left(\left[1 - \frac{l}{1-R} \right] \vee 0, F(t); -\sqrt{1 - \nu}\right) -
#'       C_2\left(\left[1 - \frac{u}{1-R} \right] \vee 0, F(t); -\sqrt{1 - \nu}\right) \right)
#' }
#' with \eqn{g(x) = \{[(1 - R) x - l] \vee (u-l)\} \wedge 1}, \eqn{C_2} being
#' the bivariate Gaussian copula, and \eqn{F} being the Exponential distribution
#' function for \eqn{\lambda}.
#'
#' @examples
#' parm <- ExtGaussian2FParam(5L, 8e-2, rho = 4e-1)
#' expected_cdo_equation(
#'   parm, times = seq(25e-2, 5, by = 25e-2), discount_factors = rep(1, 20L),
#'   recovery_rate = 0.4, lower = c(0, 0.1, 0.2, 0.35), upper = c(0.1, 0.2, 0.35, 1),
#'   coupon = c(rep(5e-2, 3L), 0), upfront = c(8e-1, 5e-1, 1e-1, -5e-2))
#' expected_cdo_equation(
#'   parm, times = seq(25e-2, 5, by = 25e-2), discount_factors = rep(1, 20L),
#'   recovery_rate = 0.4, lower = c(0, 0.1, 0.2, 0.35), upper = c(0.1, 0.2, 0.35, 1),
#'   coupon = c(rep(5e-2, 3L), 0), upfront = c(8e-1, 5e-1, 1e-1, -5e-2),
#'   method = "prob")
#' expected_cdo_equation(
#'   parm, times = seq(25e-2, 5, by = 25e-2), discount_factors = rep(1, 20L),
#'   recovery_rate = 0.4, lower = c(0, 0.1, 0.2, 0.35), upper = c(0.1, 0.2, 0.35, 1),
#'   coupon = c(rep(5e-2, 3L), 0), upfront = c(8e-1, 5e-1, 1e-1, -5e-2),
#'   method = "mc", n_sim = 1e1)
#'
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats pexp qnorm
#' @importFrom checkmate assert check_numeric
#'
#' @export
setMethod("expected_cdo_equation", "ExtGaussian2FParam",
  function(
      object, times, discount_factors, recovery_rate, lower, upper, coupon, upfront, ...,
      method = c("default", "prob", "mc")) {
    method <- match.arg(method)
    if ("default" == method) {
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

      corr <- p2P(-sqrt(1 - getNu(object)), 2L)
      x <- matrix(nrow = length(times), ncol = p)
      for (j in seq_len(ncol(x))) {
        x[, j] <- pexp(times, rate = getLambda(object)) %>%
          map_dbl(~{
            u_left <- c(pmax(1 - lower[j] / (1 - recovery_rate[j]), 0), .)
            u_right <- c(pmax(1 - upper[j] / (1 - recovery_rate[j]), 0), .)
            (1 - recovery_rate[j]) *
              (pmvnorm(rep(-Inf, 2L), qnorm(u_left), corr = corr) -
                pmvnorm(rep(-Inf, 2L), qnorm(u_right), corr = corr))
          })
      }
      out <- Rcpp__lagg_ev_cdo(
        x, times, discount_factors, recovery_rate, lower, upper, coupon, upfront)
    } else {
      out <- callNextMethod(
        object, times, discount_factors, recovery_rate, lower, upper, coupon, upfront, ...,
        method = method)
    }

    out
  })


#' @describeIn ExtGaussian2FParam-class Display the object.
#' @aliases show,ExtGaussian2FParam-method
#'
#' @inheritParams methods::show
#'
#' @export
setMethod("show", "ExtGaussian2FParam",
 function(object) {
   cat(sprintf("An object of class %s\n", classLabel(class(object))))
   if (isTRUE(validObject(object, test = TRUE))) {
     cat(sprintf("Dimension: %i\n", getDimension(object)))
     cat("Parameter:\n")
     cat(sprintf("* %s: %s\n", "Lambda", format(getLambda(object))))
     cat(sprintf("* %s: %s\n", "Rho", format(getRho(object))))
     cat(sprintf("* %s: %s\n", "Tau", format(getTau(object))))
     cat("Internal parameter:\n")
     cat(sprintf("* %s: %s\n", "Nu", format(getNu(object))))
   } else {
     cat("\t (invalid or not initialized)\n")
   }

   invisible(NULL)
  })
