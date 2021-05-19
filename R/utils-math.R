#' @keywords internal
dot <- function(x, y) {
  drop(t(x) %*% y)
}
