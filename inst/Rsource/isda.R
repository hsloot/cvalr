first <- function(x) {
  x[[1]]
}
last <- function(x) {
  x[[length(x)]]
}
adjacent_difference <- function(x) {
  c(first(x), diff(x))
}

swap_rate <- function(p, t) {
  dt <- adjacent_difference(t)

  (1 - last(p)) / (sum(dt * p))
}

zero_bond_price <- function(r, t, method = c("continuous", "discrete", "simple"), m = 1) {
  method <- match.arg(method)
  if (method == "continuous"){
    return( exp(-r * t) )
  } else if (method == "discrete") {
    stopifnot(m %% 1 == 0, 1 <= m)
    return( (1 + r / m)^(t * m) )
  } else {
    return( 1 / (1 + r * t) )
  }
}

yield <- function(p, t, method = c("continuous", "discrete", "simple"), m = 1) {
  method <- match.arg(method)
  if (method == "continuous"){
    return(-log(p) / t)
  } else if (method == "discrete") {
    stopifnot(m %% 1 == 0, 1 <= m)
    return( m * (p^(1/(t * m)) - 1) )
  } else {
   return( (1 / p - 1) / t )
  }
}

#' Calculate discount factors according to the ISDA standard model
#'
#' @param values Deposit or swap rates
#' @param maturities Maturities of the instruments
#' @param maturities.out Desired maturities the returned
#'   discount factors
#'
#' @details
#' Assumes that daycount conventions are considered in the provided
#' maturities. Also assumes that maturities less or equal to one
#' correspond to deposit rates and maturities greater than one
#' correspond to swap rates.
#' Deposit rates are assumed to be simple and directly converted
#' to discount factors; swap rates are used in an iterative approach
#' to bootstrap the  algorithm).
#'
#' __WARNING__: Might only be approximately correct; use with care.
#'
#' @returns
#' A vector with discount factors to the corresponding maturities.
discount_factors <- function(
    values, maturities, maturities.out = maturities) {
  is_deposit <- (maturities <= 1)
  discounts <- zero_bond_price(
    values[is_deposit], maturities[is_deposit], method = "simple")
  for (i in seq_along(values)[!is_deposit]) {
    brent_res <- uniroot(
      function(x) {
        tin <- maturities[1:i]
        yin <- yield(c(discounts, x), tin, method = "continuous")
        t <- seq(1/4, maturities[i], by = 1/4)
        p <- zero_bond_price(
          approx(tin, yin, xout = t)$y , t, method = "continuous")

        swap_rate(p, t) - values[i]
      },
      interval = exp(-c(1, -1) * maturities[i])
    )
    discounts <- c(discounts, brent_res$root)
  }

  discounts
}
