.onLoad <- function(libname, pkgname) { # nolint
  op <- options()
  op_cvalr <- list(
    cvalr.enable_warnings = TRUE,
    cvalr.enable_messages = TRUE
  )
  to_set <- !(names(op_cvalr) %in% names(op))
  if (any(to_set)) options(op_cvalr[to_set])
}
