#' @keywords internal
dot <- function(x, y, simplify = TRUE) {
  out <- t(x) %*% y
  if (simplify) out <- drop(out)

  out
}
