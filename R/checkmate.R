# nocov start

#' @importFrom waldo compare
checkEqual <- function(x, y, ..., tolerance = .Machine$double.eps^0.5) { # nolint
  out <- waldo::compare(x, y, ..., tolerance = tolerance)
  if (isTRUE(0L < length(out)))
    return(paste(c("x and y are not equal within tolerance", out), collapse = "\\n"))

  invisible(TRUE)
}

check_equal <- checkEqual
assertEqual <- checkmate::makeAssertionFunction(checkEqual) # nolint
assert_equal <- assertEqual

#' @importFrom checkmate check_matrix makeAssertionFunction
#' @include RcppExports.R
checkQMatrix <- function(x, ...) { # nolint
  out <- check_matrix(
      x, mode = "numeric", any.missing = FALSE, all.missing = FALSE,
      ...)
  if (!isTRUE(out))
    return(out)
  if (!isTRUE(is_qmatrix(x, tol = .Machine$double.eps^0.5)))
    return("Must be upper triangular Markov intensity matrix")

  invisible(TRUE)
}

check_qmatrix <- checkQMatrix
assertQMatrix <- checkmate::makeAssertionFunction(checkQMatrix) # nolint
assert_qmatrix <- assertQMatrix

# nocov end