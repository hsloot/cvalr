test__pcds_fn <- function(x, recovery_rate) {
  qassert(recovery_rate, "N1")

  (1 - recovery_rate) * x
}

test__cdo_fn <- function(x, recovery_rate, lower, upper) {
  qassert(recovery_rate, "N1")
  qassert(lower, "N1")
  qassert(upper, "N1")

  pmin(pmax((1 - recovery_rate) * x - lower, 0), upper - lower)
}

test__expected_pcds_equation__default <- function( # nolint
    object, times, discount_factors, recovery_rate, coupon, upfront) {
  assert_true(hasMethod("getLambda", class(object)))

  lambda <- getLambda(object)
  p <- vctrs::vec_size_common(recovery_rate, coupon, upfront)
  recovery_rate <- vctrs::vec_recycle(recovery_rate, p)
  coupon <- vctrs::vec_recycle(coupon, p)
  upfront <- vctrs::vec_recycle(upfront, p)

  x <- outer(times, recovery_rate, function(times, recovery_rate) {
    (1 - recovery_rate) * pexp(times, rate = lambda)
  })
  out <- numeric(p)
  for (i in seq_along(out)) {
    out[i] <- Rcpp__pcds_edtl(
      x[, i], times, discount_factors, recovery_rate[i], coupon[i], upfront[i])
  }

  out
}

test__expected_pcds_loss__prob <- function(object, times, recovery_rate) {
  assert_true(hasMethod("probability_distribution", class(object)))

  d <- getDimension(object)
  p <- length(recovery_rate)
  probs <- probability_distribution(object, times)
  out <- matrix(nrow = length(times), ncol = p)
  for (i in seq_len(ncol(out))) {
    vals <- sapply((0:d) / d, test__pcds_fn, recovery_rate = recovery_rate[i])
    out[, i] <- drop(t(vals) %*% probs)
  }

  drop(out)
}

test__expected_pcds_equation__prob <- function( # nolint
  object, times, discount_factors, recovery_rate, coupon, upfront) {
  assert_true(hasMethod("probability_distribution", class(object)))

  p <- vctrs::vec_size_common(recovery_rate, coupon, upfront)
  recovery_rate <- vctrs::vec_recycle(recovery_rate, p)
  coupon <- vctrs::vec_recycle(coupon, p)
  upfront <- vctrs::vec_recycle(upfront, p)

  x <- matrix(
    test__expected_pcds_loss__prob(object, times, recovery_rate),
    nrow = length(times), ncol = p)
  out <- numeric(p)
  for (i in seq_along(out)) {
    out[i] <- Rcpp__pcds_edtl(
      x[, i], times, discount_factors, recovery_rate[i], coupon[i], upfront[i])
  }

  out
}

test__expected_pcds_equation__mc <- function( # nolint
    object, times, discount_factors, recovery_rate, coupon, upfront, n_sim = 1e4L) {
  p <- vctrs::vec_size_common(recovery_rate, coupon, upfront)
  recovery_rate <- vctrs::vec_recycle(recovery_rate, p)
  coupon <- vctrs::vec_recycle(coupon, p)
  upfront <- vctrs::vec_recycle(upfront, p)
  x <- simulate_adcp(object, times, n_sim = n_sim)

  y <- matrix(nrow = n_sim, ncol = p)
  for (i in seq_len(nrow(y))) {
    for (j in seq_len(ncol(y))) {
      y[i, j] <- Rcpp__pcds_edtl(
        test__pcds_fn(x[i, ], recovery_rate[j]), times, discount_factors,
        recovery_rate = recovery_rate[j], coupon = coupon[j], upfront = upfront[j])
    }
  }

  drop(apply(y, 2L, mean))
}

test__expected_cdo_equation__gaussian <- function( # nolint
    object, times, discount_factors, recovery_rate, lower, upper, coupon, upfront) {
  assert_true(is(object, "ExtGaussian2FParam"))

  p <- vctrs::vec_size_common(recovery_rate, lower, upper, coupon, upfront)
  recovery_rate <- vctrs::vec_recycle(recovery_rate, p)
  lower <- vctrs::vec_recycle(lower, p)
  upper <- vctrs::vec_recycle(upper, p)
  coupon <- vctrs::vec_recycle(coupon, p)
  upfront <- vctrs::vec_recycle(upfront, p)

  lambda <- getLambda(object)
  nu <- getNu(object)

  corr <- copula::p2P(-sqrt(1 - nu), 2L)
  x <- matrix(nrow = length(times), ncol = p)
  for (i in seq_len(nrow(x))) {
    for (j in seq_len(ncol(x))) {
      u_left <- c(1 - pmin(lower[j] / (1 - recovery_rate[j]), 1), pexp(times[i], rate = lambda))
      u_right <- c(1 - pmin(upper[j] / (1 - recovery_rate[j]), 1), pexp(times[i], rate = lambda))

      x_left <- mvtnorm::pmvnorm(lower = rep.int(-Inf, 2L), upper = qnorm(u_left), corr = corr)
      x_right <- mvtnorm::pmvnorm(lower = rep.int(-Inf, 2L), upper = qnorm(u_right), corr = corr)

      x[i, j] <- (1 - recovery_rate[j]) * (x_left - x_right)
    }
  }
  out <- numeric(p)
  for (j in seq_along(out)) {
    out[j] <- Rcpp__cdo_edtl(
      x[, j], times, discount_factors, recovery_rate[j], lower[j], upper[j], coupon[j], upfront[j])
  }

  out
}

test__expected_cdo_loss__prob <- function(object, times, recovery_rate, lower, upper) {
  assert_true(hasMethod("probability_distribution", class(object)))

  d <- getDimension(object)
  p <- vctrs::vec_size_common(recovery_rate, lower, upper)
  recovery_rate <- vctrs::vec_recycle(recovery_rate, p)
  lower <- vctrs::vec_recycle(lower, p)
  upper <- vctrs::vec_recycle(upper, p)

  probs <- probability_distribution(object, times)
  out <- matrix(nrow = length(times), ncol = p)
  for (i in seq_len(ncol(out))) {
    vals <- sapply((0:d) / d, test__cdo_fn,
                 recovery_rate = recovery_rate[i], lower = lower[i], upper = upper[i])
    out[, i] <- drop(t(vals) %*% probs)
  }

  drop(out)
}

test__expected_cdo_equation__prob <- function( # nolint
  object, times, discount_factors, recovery_rate, lower, upper, coupon, upfront) {
  assert_true(hasMethod("probability_distribution", class(object)))

  p <- vctrs::vec_size_common(recovery_rate, lower, upper, coupon, upfront)
  recovery_rate <- vctrs::vec_recycle(recovery_rate, p)
  lower <- vctrs::vec_recycle(lower, p)
  upper <- vctrs::vec_recycle(upper, p)
  coupon <- vctrs::vec_recycle(coupon, p)
  upfront <- vctrs::vec_recycle(upfront, p)

  x <- matrix(
    test__expected_cdo_loss__prob(object, times, recovery_rate, lower, upper),
    nrow = length(times), ncol = p)
  out <- numeric(p)
  for (i in seq_along(out)) {
    out[i] <- Rcpp__cdo_edtl(
      x[, i], times, discount_factors, recovery_rate[i], lower[i], upper[i], coupon[i], upfront[i])
  }

  out
}

test__expected_cdo_equation__mc <- function( # nolint
    object, times, discount_factors, recovery_rate, lower, upper, coupon, upfront, n_sim = 1e4L) {
  p <- vctrs::vec_size_common(recovery_rate, lower, upper, coupon, upfront)
  recovery_rate <- vctrs::vec_recycle(recovery_rate, p)
  lower <- vctrs::vec_recycle(lower, p)
  upper <- vctrs::vec_recycle(upper, p)
  coupon <- vctrs::vec_recycle(coupon, p)
  upfront <- vctrs::vec_recycle(upfront, p)
  x <- simulate_adcp(object, times, n_sim = n_sim)

  y <- matrix(nrow = n_sim, ncol = p)
  for (i in seq_len(nrow(y))) {
    for (j in seq_len(ncol(y))) {
      y[i, j] <- Rcpp__cdo_edtl(
        test__cdo_fn(x[i, ], recovery_rate[j], lower[j], upper[j]), times, discount_factors,
        recovery_rate[j], lower = lower[j], upper = upper[j],
        coupon = coupon[j], upfront = upfront[j])
    }
  }

  drop(apply(y, 2L, mean))
}
