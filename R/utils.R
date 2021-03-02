#' @keywords internal
simplify2vector <- function(x) {
  if (isTRUE(is.matrix(x)) && isTRUE(1L == nrow(x) || 1L == ncol(x))) {
    x <- as.vector(x)
  }

  x
}
