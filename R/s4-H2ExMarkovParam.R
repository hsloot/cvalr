#' @include s4-H2ExCalibrationParam.R checkmate.R
NULL

# nolint start
ERR_MSG_FRACTION <- "`fraction` must be in [0,1]"
ERR_MSG_MODELS <- "`models` has wrong type"
# nolint end

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
#' @importFrom purrr map_lgl map2_lgl map_int walk
#' @importFrom checkmate qtest test_choice test_list
setValidity("H2ExMarkovParam", # nolint
  function(object) {
    if (!qtest(object@fraction, "N1(0,1)")) {
      return(ERR_MSG_FRACTION)
    }
    if (!(
        test_list(object@models, types = getModelName(object), any.missing = FALSE,
          len = length(getComposition(object)) + 1L) &&
        all(map_lgl(object@models, ~isTRUE(validObject(., test = TRUE)))) &&
        test_choice(getDimension(getGlobalModel(object)), getDimension(object)) &&
        test_choice(length(getPartitionModels(object)), length(getComposition(object))) &&
        all(map2_lgl(getPartitionModels(object), getComposition(object), ~{
          getDimension(.x) == .y
        })))) {
      return(ERR_MSG_MODELS)
    }

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
#' H2ExMarkovParam()
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
#' @param n_sim Number of samples.
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
  function(object, ..., n_sim = 1e4L) {
    fraction <- getFraction(object)
    dt_global <- simulate_dt(getGlobalModel(object), ..., n_sim = n_sim)
    dt_partition <- map(getPartitionModels(object), simulate_dt, ..., n_sim = n_sim) %>%
      reduce(cbind) %>%
      `dimnames<-`(NULL)

    pmin(1 / fraction * dt_global, 1 / (1 - fraction) * dt_partition)
  })


#' @describeIn H2ExMarkovParam-class
#'   calculates the *probability vector* for the *average default count process*
#'   and returns a matrix `x` with `dim(x) == c(getDimension(object)+1L, length(times))`.
#' @aliases probability_distribution,H2ExMarkovParam-method
#'
#' @inheritParams probability_distribution
#'
#' @examples
#' probability_distribution(ArmageddonH2ExtMO3FParam(
#'   composition = c(2L, 4L, 2L), lambda = 0.05, rho = c(3e-1, 6e-1)), 0.3)
#' probability_distribution(AlphaStableH2ExtMO3FParam(
#'   composition = c(2L, 4L, 2L), lambda = 0.05, rho = c(3e-1, 6e-1)), 0.3)
#' probability_distribution(PoissonH2ExtMO3FParam(
#'   composition = c(2L, 4L, 2L), lambda = 0.05, rho = c(3e-1, 6e-1)), 0.3)
#' probability_distribution(ExponentialH2ExtMO3FParam(
#'   composition = c(2L, 4L, 2L), lambda = 0.05, rho = c(3e-1, 6e-1)), 0.3)
#'
#' @section Probability distribution:
#' The probability of \eqn{j > i} portfolio items being defaulted at time
#' \eqn{t > s} conditioned on \eqn{i} portfolio items being defaulted at time
#' \eqn{s} is
#' \deqn{
#'   \mathbb{P}(Z_t = j \mid Z_s = i)
#'     = \sum_{i + i_0 + i_1 + \cdots + i_J = j}
#'       \delta_{i}^\top \operatorname{e}^{(t-s) Q^{(1)}} \delta_{i_1} \cdot
#'         \delta_{i_1}^\top \operatorname{e}^{(t-s) Q^{(2)}} \delta_{i_2} \cdot \cdots \cdot
#'         \delta_{i_{J-1}}^\top \operatorname{e}^{(t-s) Q^{(J)}} \delta_{i_J} \cdot
#'         \delta_{i_J}^\top \operatorname{e}^{(t-s) Q^{(0)}} \delta_{i_0} ,
#' }
#' where \eqn{Q^{(i)}, \ldots, Q^{(J)}} are the Markovian generator matrices of the partition models
#' and \eqn{Q^{(0)}} is the Markovian generator matrix of the global model.
#'
#' @importFrom stats convolve
#' @importFrom expm expm
#' @importFrom checkmate qassert
#' @importFrom purrr map map2 array_branch reduce
#' @include utils.R
#'
#' @export
setMethod("probability_distribution", "H2ExMarkovParam",
  function(object, times, ...) {
    qassert(times, "N+[0,)")
    frac <- getFraction(object)
    exq0 <- getExQMatrix(getGlobalModel(object))

    map(times, ~{
      .t <- .
      map(getPartitionModels(object), ~{
        expm(.t * (1 - frac) * getExQMatrix(.))[1L, , drop = TRUE]
      }) %>%
      reduce(~{
        pmax(convolve(.x, rev(.y), type = "open"), 0)
      }) %>%
      { as.vector(t(.) %*% expm(.t * frac * exq0)) } # nolint
    }) %>%
    map(matrix, ncol = 1L) %>%
    reduce(cbind) %>%
    `dimnames<-`(NULL)
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
    if (isTRUE(validObject(object, test = TRUE))) {
      cat(sprintf("Composition: %s = %s\n", getDimension(object),
        paste(getComposition(object), collapse = " + ")))
      cat(sprintf("Fraction: %s\n", format(getFraction(object))))
      cat("Models:\n")
      cat("* Global model\n")
      writeLines(
        paste0("\t", capture.output(show(as(getGlobalModel(object), getModelName(object))))))
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
    } else {
      cat("\t (invalid or not initialized)\n")
    }

    invisible(NULL)
    })
