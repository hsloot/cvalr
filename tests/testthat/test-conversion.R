set.seed(1623)

n <- 5e1L
d <- 75L
lambda <- 1e-1
times <- seq(25e-2, 5, by = 25e-2)

dt_mat <- matrix(rexp(n * d, rate = lambda), nrow = n, ncol = d)

R__dt2adcp <- function(x, times) { # nolint
  purrr::array_branch(x, 1L) %>%
    purrr::map(~{
      dt <- .
      map_dbl(times, ~{
        mean(dt <= .)
      }) %>%
      matrix(nrow = 1L)
    }) %>%
    reduce(rbind) %>%
    `dimnames<-`(NULL)
}

test_that("Conversion \"dt -> adcp\" works as expected", {
  expect_equal(Rcpp__dt2adcp(dt_mat, times[[1]]), R__dt2adcp(dt_mat, times[[1]]))
  expect_equal(Rcpp__dt2adcp(dt_mat, times), R__dt2adcp(dt_mat, times))
})

adcp_mat <- R__dt2adcp(dt_mat, times)
recovery_rate <- 4e-1

R__adcp2l_pcds <- function(x, recovery_rate) { # nolint
  x * (1 - recovery_rate)
}

R__adcp2l_cdo <- function(x, recovery_rate, lower, upper) { # nolint
  pmin(pmax(x * (1 - recovery_rate) - lower, 0), upper - lower)
}
discount_factors <- rep(1, length(times))

test_that("Conversion \"adcp -> peqpv (pcds)\" works as expected", {
  R__adcp2peqpv_pcds <- function(x, times, discount_factors, recovery_rate, coupon, upfront) { # nolint
    n <- nrow(x)
    p <- vctrs::vec_size_common(recovery_rate, coupon, upfront)
    recovery_rate <- vctrs::vec_recycle(recovery_rate, p)
    coupon <- vctrs::vec_recycle(coupon, p)
    upfront <- vctrs::vec_recycle(upfront, p)

    out <- matrix(nrow = n, ncol = p)
    for (j in 1:p) {
      y <- R__adcp2l_pcds(x, recovery_rate[j])
      for (k in 1:n) {
        out[k, j] <- Rcpp__pcds_edtl(
          y[k, ], times, discount_factors, recovery_rate[j], coupon[j], upfront[j])
      }
    }

    out
  }
  recovery_rate <- rep(4-1, 4L)
  coupon <- rep(9e-1, 4L)
  upfront <- rep(0, 4L)

  expect_equal(
    Rcpp__adcp2peqpv_pcds(adcp_mat, times, discount_factors, recovery_rate, coupon, upfront),
    R__adcp2peqpv_pcds(adcp_mat, times, discount_factors, recovery_rate, coupon, upfront)
  )
})

test_that("Conversion \"adcp -> peqpv (cdo)\" works as expected", {
  R__adcp2peqpv_cdo <- function( # nolint
      x, times, discount_factors, recovery_rate, lower, upper, coupon, upfront) {
    n <- nrow(x)
    p <- vctrs::vec_size_common(recovery_rate, coupon, upfront, lower, upper)
    recovery_rate <- vctrs::vec_recycle(recovery_rate, p)
    lower <- vctrs::vec_recycle(lower, p)
    upper <- vctrs::vec_recycle(upper, p)
    coupon <- vctrs::vec_recycle(coupon, p)
    upfront <- vctrs::vec_recycle(upfront, p)

    out <- matrix(nrow = n, ncol = p)
    for (j in 1:p) {
      y <- R__adcp2l_cdo(x, recovery_rate[j], lower[j], upper[j])
      for (k in 1:n) {
        out[k, j] <- Rcpp__cdo_edtl(
          y[k, ], times, discount_factors, recovery_rate[j], lower[j], upper[j], coupon[j], upfront[j])
      }
    }

    out
  }
  recovery_rate <- rep(4-1, 4L)
  coupon <- c(rep(5e-2, 3), 0)
  upfront <- c(8e-1, 5e-1, 1e-1, 0)
  lower <- c(0, 0.1, 0.2, 0.35)
  upper <- c(0.1, 0.2, 0.35, 1)

  expect_equal(
    Rcpp__adcp2peqpv_cdo(
      adcp_mat, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
    R__adcp2peqpv_cdo(
      adcp_mat, times, discount_factors, recovery_rate, lower, upper, coupon, upfront)
  )
})

test_that("Conversion \"adcp -> epd\" works as expected", {
  R__adcp2epd <- function(x, d) { # nolint
    out <- matrix(nrow = d+1, ncol = ncol(x))
    for (i in 0:d) {
      for (j in seq_len(ncol(x))) {
        out[i+1, j] <- mean(x[, j] == (i / d))
      }
    }

    out
  }

  expect_equal(Rcpp__adcp2epd(adcp_mat, d), R__adcp2epd(adcp_mat, d))

})
