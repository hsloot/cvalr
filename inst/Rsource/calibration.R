## ---- Portfolio CDS calibration ----

#' @importFrom stats pexp uniroot
portfolio_cds_calibration <- function(
    coupon, upfront,
    times, discount_factors, recovery_rate,
    expected_losses_fn = function(times, recovery_rate, lambda, ...) {
      pexp(times, rate = lambda) * (1 - recovery_rate)
    },
    ...,
    tolerance = .Machine$double.eps^0.5,
    opt_interval = -log(c(1 - tolerance, tolerance))) {
  opt_fn <- Vectorize(
    FUN = function(lambda) {
      expected_losses <- expected_losses_fn(times, recovery_rate, lambda, ...)
      portfolio_cds_equation(
        expected_losses = expected_losses,
        times = times,
        discount_factors = discount_factors,
        recovery_rate = recovery_rate,
        coupon = coupon,
        upfront = upfront
      )
    }, vectorize.args = "lambda"
  )
  if (prod(opt_fn(opt_interval)) >= 0) {
    opt_root <- opt_interval[order(abs(opt_fn(opt_interval)))[[1]]]
    opt_error <- opt_fn(opt_root)
  } else {
    opt_res <- uniroot(
      f = opt_fn, interval = opt_interval, tol = tolerance
    )
    opt_root <- opt_res$root
    opt_error <- opt_res$f.root
  }

  data.frame(
    Coupon = coupon,
    Lambda = opt_res$root, Error = opt_res$f.root
  )
}


## ---- CDO tranche expected loss ----

#' @importFrom stats pexp qnorm
#' @importFrom mvtnorm pmvnorm
cdo_tranche_expected_loss_gaussian_proxy <- function(
    times, lower, upper, recovery_rate, lambda, nu) {
  corr <- matrix(c(1, rep(-sqrt(1 - nu), 2), 1), 2, 2)
  call_expectation <- Vectorize(
    FUN = function(t, x, recovery_rate, lambda, corr) {
      if (x >= 1 - recovery_rate) {
        return(0)
      } else if (x == 0) {
        return(pexp(t, rate = lambda))
      } else {
        return(
          mvtnorm::pmvnorm(
            lower = rep(-Inf, 2),
            upper = c(
              qnorm(pmax(1 - x / (1 - recovery_rate), 0)),
              qnorm(pexp(t, rate = lambda))
            ),
            corr = corr
          )
        )
      }
    }, vectorize.args = "t"
  )

  (1 - recovery_rate) * (
    call_expectation(times, lower, recovery_rate, lambda, corr) -
    call_expectation(times, upper, recovery_rate, lambda, corr)
  )
}


#' @importFrom rmo rexmo_markovian ex_intensities_alpha_stable
cdo_tranche_expected_loss_exmo_stable_mc <- function(
    times, lower, upper, recovery_rate, lambda, nu, d = 75, n = 1e3) {
  alpha <- log2(2 - nu)
  tau <- rexmo_markovian(n, d, ex_intensities_alpha_stable(d, alpha)) / lambda
  sapply(times, function(t) {
    mean(pmin(
      pmax(
        (1 - recovery_rate) * apply(tau <= t, 1, mean) - lower,
        0
      ),
      upper - lower
    ))
  })
}


## ---- CDO tranche calibration ----

cdo_tranche_calibration <- function(
    coupon, upfront, lambda,
    times, discount_factors, lower, upper, recovery_rate,
    expected_losses_fn = function(times, lower, upper, recovery_rate, lambda, nu, ...) {
      cdo_tranche_expected_loss_gaussian_proxy(
        times, lower, upper, recovery_rate, lambda, nu
      )
    },
    ...,
    tolerance = .Machine$double.eps^0.5,
    opt_interval = c(tolerance, 1 - tolerance)) {
  opt_fn <- Vectorize(
    FUN = function(nu) {
      expected_losses <- expected_losses_fn(
        times = times, lower = lower, upper = upper,
        recovery_rate = recovery_rate, lambda = lambda, nu = nu, ...
      )
      cdo_equation(
          expected_losses = expected_losses,
          times = times,
          discount_factors = discount_factors,
          lower = lower,
          upper = upper,
          coupon = coupon,
          upfront = upfront
      )
    }, vectorize.args = "nu"
  )
  if (upper - lower == 1) {
    opt_root <- max(opt_interval)
    opt_error <- (opt_fn(min(opt_interval)) + opt_fn(max(opt_interval)))/2
  } else if (!isTRUE(prod(sign(opt_fn(opt_interval))) < 0)) {
    opt_root <- opt_interval[order(abs(opt_fn(opt_interval)))[[1]]]
    opt_error <- opt_fn(opt_root)
  } else {
    opt_res <- uniroot(
      f = opt_fn, interval = opt_interval, tol = tolerance
    )
    opt_root <- opt_res$root
    opt_error <- opt_res$f.root
  }

  data.frame(
    Coupon = coupon, Upfront = upfront, Nu = opt_root, Error = opt_error
  )
}


## ---- CDO (combined) calibration ----

cdo_combined_calibration <- function(
    coupon, upfront, lambda,
    times, discount_factors, lower, upper, recovery_rate,
    norm = function(x) {
      sum(abs(x))
    },
    expected_losses_fn = function(times, lower, upper, recovery_rate, lambda, nu, ...) {
      cdo_tranche_expected_loss_gaussian_proxy(
        times, lower, upper, recovery_rate, lambda, nu
      )
    },
    ...,
    tolerance = .Machine$double.eps^0.5,
    opt_interval = c(tolerance, 1 - tolerance)) {
  opt_fn_ <- Vectorize(
    FUN = function(nu, times, discount_factors, lower, upper, recovery_rate,
        coupon, upfront, lambda, ...) {
      expected_losses <- expected_losses_fn(
        times = times, lower = lower, upper = upper,
        recovery_rate = recovery_rate, lambda = lambda, nu = nu, ...
      )
      cdo_equation(
          expected_losses = expected_losses,
          times = times,
          discount_factors = discount_factors,
          lower = lower,
          upper = upper,
          coupon = coupon,
          upfront = upfront
      )
    },
    vectorize.args = c("lower", "upper", "recovery_rate", "coupon", "upfront", "lambda")
  )
  opt_fn <- Vectorize(
    FUN = function(nu) {
      norm(opt_fn_(nu, times, discount_factors, lower, upper, recovery_rate,
          coupon, upfront, lambda, ...))
    },
    vectorize.args = "nu"
  )
  opt_res <- optim(
    par = 0.5,
    fn = opt_fn,
    lower = min(opt_interval),
    upper = max(opt_interval),
    control = list(
      factr = tolerance
    ),
    method = "L-BFGS-B"
  )

  opt_root <- opt_res$par
  opt_error <- opt_fn_(opt_root, times, discount_factors, lower, upper,
    recovery_rate, coupon, upfront, lambda, ...)

  data.frame(
    Coupon = coupon, Upfront = upfront, Nu = opt_root, Error = opt_error
  )
}
