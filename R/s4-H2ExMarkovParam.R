#' @include s4-H2ExCalibrationParam.R checkmate.R
NULL

#' H2-exchangeable Markovian calibration parameter
#'
#' [CalibrationParam-class] for the H2-exchangeable Markovian *(average) default counting process*
#' model. Extends [H2ExCalibrationParam-class] and related to [ExMarkovParam-class].
#'
#' @slot models A list with the global and component models (of the type
#'   `ExMarkovParam-class`).
#' @slot fraction The proportion associated with the global model, see details.
#'
#' @details
#' The model is defined by the assumption that the vector of default times is defined as the
#' component-wise minimum of two vectors of the same length. The first vector is simulated from a
#' (scaled) global [ExMarkovParam-class] model and the second vector is the (scaled) conjunction of
#' independent [ExMarkovParam-class] models. The inverse scaling factors are a convex combination.
#'
#' @export H2ExMarkovParam
H2ExMarkovParam <- setClass("H2ExMarkovParam", # nolint
  contains = "H2ExCalibrationParam",
  slots = c(fraction = "numeric", models = "list"))


setGeneric("getFraction",
  function(object) {
    standardGeneric("getFraction")
  })
setMethod("getFraction", "H2ExMarkovParam",
  function(object) {
    object@fraction
  })

setGeneric("setFraction<-",
  function(object, value) {
    standardGeneric("setFraction<-")
  })
#' @importFrom checkmate qassert
setReplaceMethod("setFraction", "H2ExMarkovParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    object@fraction <- value

    invisible(object)
  })

setGeneric("getModelName",
  function(object) {
    standardGeneric("getModelName")
  })
setMethod("getModelName", "H2ExMarkovParam",
  function(object) {
    "ExMarkovParam"
  })

setGeneric("getModels",
  function(object) {
    standardGeneric("getModels")
  })
setMethod("getModels", "H2ExMarkovParam",
  function(object) {
    object@models
  })

setGeneric("setModels<-",
  function(object, value) {
    standardGeneric("setModels<-")
  })
#' @importFrom purrr map_lgl map_int
#' @importFrom checkmate test_class assert_choice
setReplaceMethod("setModels", "H2ExMarkovParam",
  function(object, value) {
    assert_list(value, types = getModelName(object), any.missing = FALSE, min.len = 2L)
    dim <- getDimension(value[[1]])
    composition <- map_int(value[-1], getDimension)
    assert_choice(sum(composition), dim)
    setComposition(object) <- composition
    object@models <- value

    invisible(object)
  })

setGeneric("getGlobalModel",
  function(object) {
    standardGeneric("getGlobalModel")
  })
#' @importFrom checkmate assert_list
setMethod("getGlobalModel", "H2ExMarkovParam",
  function(object) {
    assert_list(object@models, min.len = 1L)
    object@models[[1L]]
  })

setGeneric("getPartitionModels",
  function(object) {
    standardGeneric("getPartitionModels")
  })
#' @importFrom checkmate assert_list
setMethod("getPartitionModels", "H2ExMarkovParam",
  function(object) {
    assert_list(object@models, min.len = 2L)
    object@models[-1L]
  })


#' @importFrom methods is
#' @importFrom purrr map_lgl map_int walk
#' @importFrom checkmate qassert assert_true assert_list
setValidity("H2ExMarkovParam",
  function(object) {
    qassert(object@fraction, "N1(0,1)")
    assert_list(
      object@models, types = getModelName(object), any.missing = FALSE,
      len = length(object@composition) + 1L)
    walk(object@models, validObject)
    assert_true(getDimension(object@models[[1]]) == object@dim)
    assert_true(all(map_int(object@models[-1], getDimension) == object@composition))

    invisible(TRUE)
  })


#' @describeIn H2ExMarkovParam-class Constructor
#' @aliases initialize,H2ExMarkovParam-method
#' @aliases initialize,H2ExMarkovParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param fraction The proportion associated with the global model, see details.
#' @param models A list with the global and component models (of the type
#'   `ExMarkovParam-class`).
#' @param ... Pass-through parameters.
#'
#' @examples
#' composition <- c(2L, 4L, 2L)
#' d <- sum(composition)
#' model_global <- ExMarkovParam(rmo::exQMatrix(rmo::AlphaStableBernsteinFunction(0.4), d))
#' model_partition <- purrr::map(composition, ~{
#'   ExMarkovParam(rmo::exQMatrix(rmo::AlphaStableBernsteinFunction(0.5), .x))
#'   })
#' models <- c(list(model_global), model_partition)
#' H2ExMarkovParam(fraction = 0.4, models = models)
#'
#' @importFrom purrr map_lgl map_int reduce
#' @importFrom checkmate assert_true assert_list
setMethod("initialize", "H2ExMarkovParam",
  function(.Object, fraction, models) { # nolint
    if (!missing(models) && !missing(fraction)) {
      setFraction(.Object) <- fraction
      setModels(.Object) <- models

      validObject(.Object)
    }

    invisible(.Object)
  })


#' @describeIn H2ExMarkovParam-class
#'    simulates the vector of *default times* and returns a matrix `x` with
#'    `dim(x) == c(n_sim, getDimension(object))`.
#' @aliases simulate_dt,H2ExMarkovParam-method
#'
#' @inheritParams simulate_dt
#'
#' @section Simulation:
#' The default times are sampled using the stochastic representation described in details.
#'
#' @examples
#' composition <- c(2L, 4L, 2L)
#' d <- sum(composition)
#' model_global <- ExMarkovParam(rmo::exQMatrix(rmo::AlphaStableBernsteinFunction(0.4), d))
#' model_partition <- purrr::map(composition, ~{
#'   ExMarkovParam(rmo::exQMatrix(rmo::AlphaStableBernsteinFunction(0.5), .x))
#'   })
#' models <- c(list(model_global), model_partition)
#' parm <- H2ExMarkovParam(fraction = 0.4, models = models)
#' simulate_dt(parm, n_sim = 5e1)
#'
#' @importFrom purrr map reduce
#'
#' @include utils.R
setMethod("simulate_dt", "H2ExMarkovParam",
  function(object, ...) {
    fraction <- getFraction(object)
    dt_global <- simulate_dt(getGlobalModel(object), ...)
    dt_partition <- map(getPartitionModels(object), simulate_dt, ...) %>%
      reduce(cbind) %>%
      `dimnames<-`(NULL)

    pmin(1 / fraction * dt_global, 1 / (1 - fraction) * dt_partition)
  })


#' @describeIn H2ExMarkovParam-class Display the object.
#' @aliases show,H2ExMarkovParam-method
#'
#' @importFrom utils capture.output
#' @importFrom purrr map compose flatten_chr
#'
#' @export
setMethod("show", "H2ExMarkovParam",
  function(object) {
    cat(sprintf("An object of class %s\n", classLabel(class(object))))
    cat(sprintf("Composition: %s = %s\n", getDimension(object),
      paste(getComposition(object), collapse = " + ")))
    cat(sprintf("Fraction: %s\n", format(getFraction(object))))
    cat("Models:\n")
    cat("* Global model\n")
    writeLines(paste0("\t", capture.output(show(as(getGlobalModel(object), getModelName(object))))))
    cat("* Partition models:\n")
    to_list_item <- function(x) {
      out <- rep("  ", length(x))
      out[[1]] <- "- "

      paste0(out, x)
    }
    getPartitionModels(object) %>%
      map(compose(to_list_item, ~capture.output(show(.)), ~as(., getModelName(object)))) %>%
      flatten_chr %>%
      paste0("\t", .) %>%
      writeLines

    invisible(NULL)
    })
