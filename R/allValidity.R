#' @include allClass.R allGeneric.R
NULL

setValidity("CalibrationParam",
  function(object) {
    stopifnot(object@dim >= 0L)

    invisible(TRUE)
  })

setValidity("ExMarkovParam",
  function(object) {
    stopifnot(object@dim+1 == nrow(object@qmatrix),
      nrow(object@qmatrix) == ncol(object@qmatrix),
      all(object@qmatrix[lower.tri(object@qmatrix)] == 0),
      all(object@qmatrix[upper.tri(object@qmatrix)] >= 0),
      isTRUE(all.equal(rep(0, nrow(object@qmatrix)), apply(object@qmatrix, 1, sum),
                tol = .Machine$double.eps^0.5)))

    invisible(TRUE)
  })

setValidity("ExMOParam",
  function(object) {
    stopifnot(object@dim == length(object@ex_intensities),
      all(object@ex_intensities >= 0),
      any(object@ex_intensities > 0))

    invisible(TRUE)
  })

setValidity("ExtMOParam",
  function(object) {
    stopifnot(is(object@bf, "BernsteinFunction"))

    invisible(TRUE)
  })

setValidity("ExtMO2FParam",
  function(object) {
    stopifnot(1L == length(object@lambda), object@lambda > 0,
      1L == length(object@nu),
      is(object@bf, "ScaledBernsteinFunction"))

    invisible(TRUE)
  })

#' @importFrom rmo valueOf ScaledBernsteinFunction SumOfBernsteinFunctions
#'   LinearBernsteinFunction ConstantBernsteinFunction
setValidity("CuadrasAugeExtMO2FParam",
  function(object) {
    stopifnot(is(object@bf, "ScaledBernsteinFunction"),
      isTRUE(all.equal(object@lambda, object@bf@scale)),
      is(object@bf@original, "SumOfBernsteinFunctions"),
      is(object@bf@original@first, "LinearBernsteinFunction"),
      is(object@bf@original@second, "ConstantBernsteinFunction"),
      isTRUE(all.equal(1, valueOf(object@bf@original, 1, 0L))),
      isTRUE(all.equal(
        object@nu, invAlpha(object, 2 - valueOf(object@bf@original, 2, 0L)))))
  })

  #' @importFrom rmo valueOf ScaledBernsteinFunction SumOfBernsteinFunctions
  #'   AlphaStableBernsteinFunction
setValidity("AlphaStableExtMO2FParam",
  function(object) {
    stopifnot(is(object@bf, "ScaledBernsteinFunction"),
      isTRUE(all.equal(object@lambda, object@bf@scale)),
      is(object@bf@original, "AlphaStableBernsteinFunction"),
      isTRUE(all.equal(1, valueOf(object@bf@original, 1, 0L))),
      isTRUE(all.equal(
        object@nu, invAlpha(object, 2 - valueOf(object@bf@original, 2, 0L)))))
  })

  #' @importFrom rmo valueOf ScaledBernsteinFunction SumOfBernsteinFunctions
  #'   LinearBernsteinFunction PoissonBernsteinFunction
setValidity("PoissonExtMO2FParam",
  function(object) {
    stopifnot(is(object@bf, "ScaledBernsteinFunction"),
      isTRUE(all.equal(object@lambda, object@bf@scale)),
      is(object@bf@original, "SumOfBernsteinFunctions"),
      is(object@bf@original@first, "LinearBernsteinFunction"),
      is(object@bf@original@second, "PoissonBernsteinFunction"),
      isTRUE(all.equal(1, valueOf(object@bf@original, 1, 0L))),
      isTRUE(all.equal(
        object@nu, invAlpha(object, 2 - valueOf(object@bf@original, 2, 0L)))))
  })

  #' @importFrom rmo valueOf ScaledBernsteinFunction SumOfBernsteinFunctions
  #'   LinearBernsteinFunction ExponentialBernsteinFunction
setValidity("ExponentialExtMO2FParam",
  function(object) {
    stopifnot(is(object@bf, "ScaledBernsteinFunction"),
      isTRUE(all.equal(object@lambda, object@bf@scale)),
      is(object@bf@original, "SumOfBernsteinFunctions"),
      is(object@bf@original@first, "LinearBernsteinFunction"),
      is(object@bf@original@second, "ExponentialBernsteinFunction"),
      isTRUE(all.equal(1, valueOf(object@bf@original, 1, 0L))),
      isTRUE(all.equal(
        object@nu, invAlpha(object, 2 - valueOf(object@bf@original, 2, 0L)))))
  })

setValidity("ExtGaussian2FParam",
  function(object) {
    stopifnot(1L == length(object@lambda), object@lambda > 0,
      1L == length(object@nu), 0 <= object@nu, object@nu <= 1)

    invisible(TRUE)
  })

setValidity("ExtArch2FParam",
  function(object) {
    stopifnot(1L == length(object@lambda), object@lambda > 0,
      1L == length(object@nu))

    invisible(TRUE)
  })
