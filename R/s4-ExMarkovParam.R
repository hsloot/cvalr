#' @include s4-CalibrationParam.R checkmate.R
NULL

#' Exchangeable Markovian calibration parameter
#'
#' Calibration parameter class for the general exchangeable model with a
#' Markovian *default counting process*.
#'
#' @slot qmatrix The \eqn{(d+1) \times (d+1)} Markov generator matrix of the
#'   default counting process.
#'
#' @details
#' The probability of \eqn{j > i} portfolio items being defaulted at time
#' \eqn{t > s} conditioned on \eqn{i} portfolio items being defaulted at time
#' \eqn{s} is
#'
#' \deqn{
#'   \mathbb{P}(Z_t = j \mid Z_s = i)
#'     = \delta_{i}^\top \operatorname{e}^{(t-s) Q} \delta_{j} .
#' }
#'
#' @export ExMarkovParam
ExMarkovParam <- setClass("ExMarkovParam", # nolint
  contains = "CalibrationParam",
  slots = c(ex_qmatrix = "matrix"))


setGeneric("getExQMatrix",
  function(object) {
    standardGeneric("getExQMatrix")
  })
setGeneric("setExQMatrix<-",
  function(object, value) {
    standardGeneric("setExQMatrix<-")
  })


setMethod("getExQMatrix", "ExMarkovParam",
  function(object) {
    object@ex_qmatrix
  })

setReplaceMethod("setExQMatrix", "ExMarkovParam",
  function(object, value) {
    assert_exqmatrix(value, min.rows = 1L, min.cols = 1L)

    dim <- nrow(value)-1
    setDimension(object) <- dim
    object@ex_qmatrix <- value

    invisible(object)
  })


setValidity("ExMarkovParam",
  function(object) {
    assert_exqmatrix(object@ex_qmatrix,
      nrows = object@dim+1L, ncols = object@dim + 1L)

    invisible(TRUE)
  })


#' @describeIn ExMarkovParam-class Constructor
#' @aliases initialize,ExMarkovParam-method
#' @aliases initialize,ExMarkovParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ex_qmatrix (Exchangeable) Q-matrix of the default counting process.
#'
#' @examples
#' ExMarkovParam(
#'  ex_qmatrix = matrix(
#'    c(-0.07647059, 0, 0, 0.05294118, -0.05, 0, 0.02352941, 0.05, 0),
#'    nrow = 3L, ncol = 3L))
setMethod("initialize", "ExMarkovParam",
  function(.Object, ex_qmatrix) { # nolint
    if (!missing(ex_qmatrix)) {
      setExQMatrix(.Object) <- ex_qmatrix
      validObject(.Object)
    }

    invisible(.Object)
  })
