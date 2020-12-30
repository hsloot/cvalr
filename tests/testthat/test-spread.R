lambda <- 0.1
rate <- 0.005
recovery_rate <- 0.4
times <- seq(0, 2, by = 0.25)
expected_default_count <- 1 - exp(-lambda * times)
discount_factors <- exp(-rate * times)

test_that("portfolio_cds_spread is calculated correctly", {
  expected_losses <- (1 - recovery_rate) * expected_default_count
  expected_nominals <- 1 - expected_losses / (1 - recovery_rate)
  expect_equal(
    portfolio_cds_spread(
      expected_losses, times, discount_factors, recovery_rate
    ),
    sum(discount_factors[-1] * diff(expected_losses)) /
      sum(
        discount_factors[-1] * diff(times) *
          0.5 * (expected_nominals[-1] + expected_nominals[-9])
      )
  )
})

lower <- 0.03
upper <- 0.1
rho <- 0.1
corr <- matrix(c(1, rep(-sqrt(1 - rho), 2), 1), 2, 2)
tmp1 <- expected_default_count
if (lower > 0) {
  tmp1 <- sapply(
    expected_default_count,
    function(x) {
      mvtnorm::pmvnorm(
        lower = rep(-Inf, 2),
        upper = c(
          -qnorm(pmin(lower/(1-recovery_rate), 1)),
          qnorm(x)
        ),
        corr=corr
      )
    }
  )
}
tmp2 <- sapply(
  expected_default_count,
  function(x) {
    mvtnorm::pmvnorm(
      lower = rep(-Inf, 2),
      upper = c(
        -qnorm(pmin(upper/(1-recovery_rate), 1)),
        qnorm(x)
      ),
      corr=corr
    )
  }
)
expected_losses <- (1 - recovery_rate) * (tmp1 - tmp2)
expected_nominals <- upper - lower - expected_losses

test_that("upfront payment is calculated correctly", {
  spread <- 0.05
  expect_equal(
    upfront_payment(
      expected_losses, times, discount_factors, lower, upper, spread
    ),
    (sum(discount_factors[-1] * diff(expected_losses)) - spread *
      sum(
        discount_factors[-1] * diff(times) *
          0.5 * (expected_nominals[-1] + expected_nominals[-9])
      )) / (upper - lower)
  )
})

test_that("upfront spread is calculated correctly", {
  expect_equal(
    upfront_spread(
      expected_losses, times, discount_factors, lower, upper
    ),
    sum(discount_factors[-1] * diff(expected_losses)) /
      sum(
        discount_factors[-1] * diff(times) *
          0.5 * (expected_nominals[-1] + expected_nominals[-9])
      )
  )
})
