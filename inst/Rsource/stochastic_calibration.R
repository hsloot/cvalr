#' @importFrom rmo rexmo_markovian ex_intensities_alpha_stable
cdo_tranche_expected_loss_exmo_alpha_stable_mc <- function(
    times, lower, upper, recovery_rate, lambda, nu, d = 75, n = 1e3) {
  alpha <- log2(2 - nu)
  tau <- rmo::rexmo_markovian(n, d, ex_intensities_alpha_stable(d, alpha)) / lambda
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

cdo_fitness_mc <- function(
    nu,
    lambda,
    times,
    discount_factors,
    lower, upper, recovery_rate,
    coupon, upfront,
    d = 75, n = 1e3,
    expected_losses_fn = cdo_tranche_expected_loss_exmo_alpha_stable_mc,
    norm_fn = function(x) { sqrt(sum(x^2)) }) {
  cdo_equations <- Vectorize(
    function(nu, lambda, lower, upper, recovery_rate, coupon, upfront) {
      expected_losses <- expected_losses_fn(times, lower, upper, recovery_rate,
        lamnda, nu, d, n)
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
    vectorize.args = c("lower", "upper", "recovery_rate", "coupon", "upfront")
  )

  -norm_fn(cdo_equations(nu, lambda, upper, recovery_rate, coupon, upfront))
}

cdo_calibration
