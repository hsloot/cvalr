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
checkExQMatrix <- function(x, ...) { # nolint
  out <- check_matrix(
      x, mode = "numeric", any.missing = FALSE, all.missing = FALSE,
      ...)
  if (!isTRUE(out))
    return(out)
  if (!isTRUE(is_exqmatrix(x, tol = .Machine$double.eps^0.5)))
    return("Must be upper triangular Markov intensity matrix")

  invisible(TRUE)
}

check_exqmatrix <- checkExQMatrix
assertExQMatrix <- checkmate::makeAssertionFunction(checkExQMatrix) # nolint
assert_exqmatrix <- assertExQMatrix

# nocov end
