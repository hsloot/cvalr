test__eddl <- function(l, df) {
  sum(df * diff(c(0, l)))
}

test__pcds_ddl <- function(l, df, recovery_rate) {
  l <- (1 - recovery_rate) * l
  test__eddl(l, df)
}

test__cdo_ddl <- function(l, df, recovery_rate, lower, upper) {
  l <- pmin(pmax((1 - recovery_rate) * l - lower, 0), upper - lower)
  test__eddl(l, df)
}

test__pcds_edpl <- function(l, times, df, recovery_rate, coupon, upfront) {
  nom <- 1 - c(0, l) / (1 - recovery_rate)
  nom <- (nom[-1] + nom[-length(nom)]) / 2
  upfront + coupon * sum(df * diff(c(0, times)) * nom)
}

test__pcds_dpl <- function(l, times, df, recovery_rate, coupon, upfront) {
  l <- (1 - recovery_rate) * l
  test__pcds_edpl(l, times, df, recovery_rate, coupon, upfront)
}

test__cdo_edpl <- function(l, times, df, recovery_rate, lower, upper, coupon, upfront) {
  nom <- (upper - lower) - c(0, l)
  nom <- (nom[-1] + nom[-length(nom)]) / 2
  upfront * (upper - lower) + coupon * sum(df * diff(c(0, times)) * nom)
}

test__cdo_dpl <- function(l, times, df, recovery_rate, lower, upper, coupon, upfront) {
  l <- pmin(pmax((1 - recovery_rate) * l - lower, 0), upper - lower)
  test__cdo_edpl(l, times, df, recovery_rate, lower, upper, coupon, upfront)
}

test__pcds_dtl <- function(l, times, df, recovery_rate, coupon, upfront) {
  test__pcds_ddl(l, df, recovery_rate) -
    test__pcds_dpl(l, times, df, recovery_rate, coupon, upfront)
}

test__cdo_dtl <- function(l, times, df, recovery_rate, lower, upper, coupon, upfront) {
  test__cdo_ddl(l, df, recovery_rate, lower, upper) -
    test__cdo_dpl(l, times, df, recovery_rate, lower, upper, coupon, upfront)
}

test__pcds_edtl <- function(l, times, df, recovery_rate, coupon, upfront) {
  test__eddl(l, df) -
    test__pcds_edpl(l, times, df, recovery_rate, coupon, upfront)
}

test__cdo_edtl <- function(l, times, df, recovery_rate, lower, upper, coupon, upfront) {
  test__eddl(l, df) -
    test__cdo_edpl(l, times, df, recovery_rate, lower, upper, coupon, upfront)
}
