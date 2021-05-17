#' @include s4-CalibrationParam.R checkmate.R
NULL

#' Virtual super-class for H2-exchangeable calibration parameters
#'
#' [H2ExCalibrationParam-class] provides a simple interface to calculate
#' expected values and pricing equations for portfolio CDS' and CDO's with
#' different hierachically-exchangeable models with 2 levels of hierarchy.
#' Extends [CalibrationParam-class].
#'
#' @slot composition Positive integerish vector for the component-composition.
#'
#' @details
#' A *H2-Exchangeable distribution* is a distribution for which a partition
#' \eqn{\mathcal{P} = \{ p_1, \ldots p_k \}} with \eqn{\dot{\bigcup}_{j} p_j = \{ 1 ,
#' \ldots, d \}} and a \eqn{\sigma}-algebra \eqn{\mathcal{G}} can be found such
#' that conditioned on \eqn{\mathcal{G}}, the subvectors corresponding the
#' elements of \eqn{\mathcal{P}} are independent and exchangeable.
#' We assume for simplicity that elements from \eqn{\mathcal{P}} contain
#' adjacent integer values such that \eqn{\mathcal{P}} can be identified with  a
#' composition \eqn{\mathcal{C} = \{ c_1, \ldots c_k \}} with \eqn{d = c_1 +
#' \cdots + c_k} and \eqn{\lvert p_i \rvert = c_i}.
#'
#' @export
setClass("H2ExCalibrationParam",
  contains = c("CalibrationParam", "VIRTUAL"),
  slots = c(composition = "integer"))


setGeneric("getComposition",
  function(object) {
    standardGeneric("getComposition")
  })
setMethod("getComposition", "H2ExCalibrationParam",
  function(object) {
    object@composition
  })

setGeneric("setComposition<-",
  function(object, value) {
    standardGeneric("setComposition<-")
  })
#' @importFrom checkmate qassert
setReplaceMethod("setComposition", "H2ExCalibrationParam",
  function(object, value) {
    qassert(value, "X+[1,)")
    qassert(sum(value), "X1[2,)")
    object@composition <- as.integer(value)
    setDimension(object) <- as.integer(sum(value))

    invisible(object)
  })

setGeneric("getPartition",
  function(object) {
    standardGeneric("getPartition")
  })
#' @importFrom utils head
#' @importFrom purrr map2
setMethod("getPartition", "H2ExCalibrationParam",
  function(object) {
    composition <- getComposition(object)
    map2(c(0, cumsum(head(composition, -1L))), composition, ~{
        .x + 1:.y
      })
  })


#' @importFrom checkmate qassert assert_choice
setValidity("H2ExCalibrationParam",
  function(object) {
    qassert(object@composition, "I+[1,)")
    assert_choice(sum(object@composition), object@dim)

    invisible(TRUE)
  })
