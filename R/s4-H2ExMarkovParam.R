#' @include s4-H2ExCalibrationParam.R checkmate.R
NULL

#' H2-exchangeable Markovian calibration parameter
#'
#' Calibration parameter class for the general 2-level
#' hierarchically-exchangeable model with Markovian *default counting processes*
#' in the exchangeable sub-components.
#'
#' @slot models A list with the global and component models (of the type
#'   `ExMarkovParam-class`)
#' @slot fraction The proportion associated with the global model, see details.
#'
#' @details
#' We assume that \eqn{\tau} has the stochastic representation to be the
#' component-wise minimum of a global exchangeable Markovian-vector
#' \eqn{\tau^{(0)}} and a vector \eqn{(\tau^{(1)}, \ldots, \tau^{(J)})} with
#' independent exchangeable Markovian-vector sub-vectors \eqn{\tau^{(j)}}.
#'
#' @export H2ExMarkovParam
H2ExMarkovParam <- setClass("H2ExMarkovParam", # nolint
  contains = "H2ExCalibrationParam",
  slots = c(models = "list", fraction = "numeric"))

setGeneric("getFraction",
  function(object) {
    standardGeneric("getFraction")
  })
setGeneric("setFraction<-",
  function(object, value) {
    standardGeneric("setFraction<-")
  })

setGeneric("getModelName",
  function(object) {
    standardGeneric("getModelName")
  })

setGeneric("getModels",
  function(object) {
    standardGeneric("getModels")
  })
setGeneric("setModels<-",
  function(object, value) {
    standardGeneric("setModels<-")
  })


setMethod("getFraction", "H2ExMarkovParam",
  function(object) {
    object@fraction
  })
#' @importFrom checkmate qassert
setReplaceMethod("setFraction", "H2ExMarkovParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    object@fraction <- value

    invisible(object)
  })

setMethod("getModels", "H2ExMarkovParam",
  function(object) {
    object@models
  })
#' @importFrom purrr map_lgl map_int map2
#' @importFrom checkmate test_class
setReplaceMethod("setModels", "H2ExMarkovParam",
  function(object, value) {
    assert_true(all(map_lgl(value, test_class, classes = getModelName(object))))
    dims <- map_int(value, getDimension)
    assert_true(dims[[1]] == sum(dims[-1]))
    partition <- map2(
      dims[-1], cumsum(c(0, dims[2:(length(dims)-1)])), ~{
        .y + 1:.x
      })
    setPartition(object) <- partition
    object@models <- value

    invisible(object)
  })


#' @importFrom methods is
#' @importFrom purrr map_lgl map2_lgl
#' @importFrom checkmate qassert assert_true
setValidity("H2ExMarkovParam",
  function(object) {
    qassert(object@fraction, "N1(0,1)")
    assert_true(all(map_lgl(object@models, ~is(.x, getModelName(object)))))
    assert_true(getDimension(object@models[[1]]) == getDimension(object))
    assert_true(length(object@models) == length(object@partition) + 1L)
    assert_true(all(map2_lgl(object@models[-1], object@partition, ~{
      getDimension(.x) == length(.y)
      })))

    invisible(TRUE)
  })


#' @importFrom purrr map_lgl map_int reduce
#' @importFrom checkmate assert_true
setMethod("initialize", "H2ExMarkovParam",
  function(.Object, models = NULL, fraction = NULL) { # nolint
    if (!missing(models) && !missing(fraction)) {
      assert_true(all(map_lgl(models, is, class = "CalibrationParam")))
      dims <- map_int(models, getDimension)
      if (length(dims) > 1L) {
        partition <- reduce(
          dims[-1], ~{
            c(.x, list(max(c(0L, unlist(.x))) + 1:.y))
          }, .init = list())
      } else {
        partition <- list(1:dims)
      }

      .Object@dim <- dims[[1]]
      .Object@partition <- partition
      .Object@models <- models
      .Object@fraction <- fraction

      validObject(.Object)
    }

    invisible(.Object)
  })


setMethod("getModelName", "H2ExMarkovParam",
  function(object) {
    "ExMarkovParam"
  })


#' @describeIn H2ExMarkovParam-class
#'    simulates the default times \eqn{(\tau_1, \ldots, \tau_d)} and returns a
#'    matrix `x` with `nrow(x) == n_sim` and `ncol(x) == dim(object)` if
#'    `dim(object) > 1L` and a vector `x` with `length(x) == n_sim` otherwise.
#' @aliases simulate_dt,H2ExMarkovParam-method
#'
#' @inheritParams simulate_dt
#'
#' @examples
#' parm <- ExponentialH2ExtMO3FParam(
#'   partition = list(1:2, 3:6, 7:8), fraction = 0.4,
#'   lambda = 8e-2, rho = c(0.2, 0.7))
#' simulate_dt(parm, n_sim = 5e1)
#'
#' @include utils.R
setMethod("simulate_dt", "H2ExMarkovParam",
  function(object, ...) {
    fraction <- object@fraction
    tmp0 <- simulate_dt(object@models[[1]], ...)
    tmp1 <- reduce(map(object@models[-1], simulate_dt, ...), cbind)
    out <- pmin(1 / fraction * tmp0, 1 / (1 - fraction) * tmp1)

    simplify2vector(out)
  })
