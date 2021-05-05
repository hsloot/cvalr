#' @keywords internal
simplify2vector <- function(x) {
  if (isTRUE(is.matrix(x)) && isTRUE(1L == nrow(x) || 1L == ncol(x))) {
    x <- as.vector(x)
  }

  x
}

#' @keywords internal
adjacent_differences <- function(x) {
  c(x[[1]], diff(x))
}
