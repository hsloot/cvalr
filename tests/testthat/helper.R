cds_expected_loss <- function(t, lambda, recovery_rate) {
  sapply(t, function(.t) (1 - recovery_rate) * pexp(.t, rate = lambda))
}

cdo_tranche_expected_loss_gaussian <- function(t, lambda, rho, lower, upper, recovery_rate) {
  sapply(t,
    function(.t, lambda, rho, lower, upper, recovery_rate) {
      corr <- matrix(c(1, rep(-sqrt(1 - rho), 2), 1), 2, 2)
      tmp1 <- pexp(.t, rate = lambda)
      if (lower > 0) {
        tmp1 <- mvtnorm::pmvnorm(
          lower = rep(-Inf, 2),
          upper = c(
            -qnorm(pmin(lower / (1 - recovery_rate), 1)),
            qnorm(pexp(.t, rate = lambda))
          ),
          corr = corr
        )
      }
      tmp2 <- mvtnorm::pmvnorm(
        lower = rep(-Inf, 2),
        upper = c(
          -qnorm(pmin(upper / (1 - recovery_rate), 1)),
          qnorm(pexp(.t, rate = lambda))
        ),
        corr = corr
      )
      (1 - recovery_rate) * (tmp1 - tmp2)
    },
    lambda = lambda, rho = rho, lower = lower, upper = upper, recovery_rate = recovery_rate
  )
}
