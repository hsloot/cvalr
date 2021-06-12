#' @include s4-H2ExMarkovParam.R
NULL

#' H2-Exchangeable Marshall--Olkin calibration parameter
#'
#' [CalibrationParam-class] for the H2-exchangeable Marshall-Olkin *(average) default counting
#' process* model. Extends [H2ExMarkovParam-class] and related to [ExMOParam-class].
#'
#' @export H2ExMOParam
H2ExMOParam <- setClass("H2ExMOParam", # nolint
  contains = "H2ExMarkovParam")

setMethod("getModelName", "H2ExMOParam",
  function(object) {
    "ExMOParam"
  })

#' @describeIn H2ExMOParam-class
#'    simulates the vector of *default times* and returns a matrix `x` with
#'    `dim(x) == c(n_sim, getDimension(object))`.
#' @aliases simulate_dt,H2ExMOParam-method
#' @param n_sim Number of samples.
#'
#' @inheritParams simulate_dt
#' @param n_sim Number of samples.
#'
#' @section Simulation:
#' The default times are sampled using the stochastic representation described in details.
#'
#' @examples
#' composition <- c(2L, 4L, 2L)
#' d <- sum(composition)
#' model_global <- ExMOParam(rmo::exIntensities(rmo::AlphaStableBernsteinFunction(0.4), d))
#' model_partition <- purrr::map(composition, ~{
#'   ExMOParam(rmo::exIntensities(rmo::AlphaStableBernsteinFunction(0.5), .x))
#'   })
#' models <- c(list(model_global), model_partition)
#' parm <- H2ExMOParam(fraction = 0.4, models = models)
#' simulate_dt(parm, n_sim = 1e1L)
#'
#' @importFrom purrr map reduce
#'
#' @include utils.R
setMethod("simulate_dt", "H2ExMOParam",
  function(object, ..., n_sim = 1e4L) {
    Rcpp__rh2exmo_markovian_dt(n_sim, getFraction(object), getModels(object))
  })

#' @describeIn H2ExMOParam-class
#'   simulates the *average default counting process* and returns a
#'   matrix `x` with `dim(x) == c(n_sim, length(times))`.
#' @aliases simulate_adcp,H2ExMOParam-methods
#'
#' @examples
#' composition <- c(2L, 4L, 2L)
#' d <- sum(composition)
#' model_global <- ExMOParam(rmo::exIntensities(rmo::AlphaStableBernsteinFunction(0.4), d))
#' model_partition <- purrr::map(composition, ~{
#'   ExMOParam(rmo::exIntensities(rmo::AlphaStableBernsteinFunction(0.5), .x))
#'   })
#' models <- c(list(model_global), model_partition)
#' parm <- H2ExMOParam(fraction = 0.4, models = models)
#' simulate_adcp(parm, 1, n_sim = 1e1L)
#' simulate_adcp(parm, seq(25e-2, 5, by = 25e-2), n_sim = 1e1L)
#'
#' @importFrom stats rexp
#' @include RcppExports.R
#'
#' @export
setMethod("simulate_adcp", "H2ExMOParam",
  function(object, times, ..., n_sim = 1e4L) {
    qassert(times, "N+[0,)")

      Rcpp__rh2exmo_markovian_adcp(n_sim, times, getFraction(object), getModels(object))
  })
