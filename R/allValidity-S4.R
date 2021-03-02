#' @include allClass-S4.R allGeneric-S4.R checkmate.R
NULL

#' @importFrom checkmate qassert
setValidity("CalibrationParam",
  function(object) {
    qassert(object@dim, "I1(0,)")

    invisible(TRUE)
  })

setValidity("ExMarkovParam",
  function(object) {
    assert_exqmatrix(object@ex_qmatrix,
      nrows = object@dim+1L, ncols = object@dim + 1L)

    invisible(TRUE)
  })

#' @importFrom checkmate qassert assert_numeric
setValidity("ExMOParam",
  function(object) {
    assert_numeric(object@ex_intensities, lower = 0, len = object@dim)
    qassert(max(object@ex_intensities), "N1(0,)")

    invisible(TRUE)
  })

#' @importFrom checkmate assert_class
setValidity("ExtMOParam",
  function(object) {
    assert_class(object@bf, "BernsteinFunction")

    invisible(TRUE)
  })

#' @importFrom rmo valueOf ScaledBernsteinFunction
#' @importFrom checkmate assert qassert check_choice check_class
setValidity("ExtMO2FParam",
  function(object) {
    qassert(object@lambda, "N1(0,)")
    qassert(object@nu, "N1")
    assert(combine = "and",
      check_class(object@bf, "ScaledBernsteinFunction"),
      check_choice(object@bf@scale, object@lambda),
      check_equal(1, valueOf(object@bf@original, 1, 0L)),
      check_equal(
        object@nu, invAlpha(object, 2 - valueOf(object@bf@original, 2, 0L))))

    invisible(TRUE)
  })

#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions
#'   LinearBernsteinFunction ConstantBernsteinFunction
#' @importFrom checkmate assert check_choice check_class
setValidity("CuadrasAugeExtMO2FParam",
  function(object) {
    assert(combine = "and",
      check_choice(object@bf@scale, object@lambda),
      check_class(object@bf@original, "SumOfBernsteinFunctions"),
      check_class(object@bf@original@first, "LinearBernsteinFunction"),
      check_class(object@bf@original@second, "ConstantBernsteinFunction"))

    invisible(TRUE)
  })

#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions
#'   AlphaStableBernsteinFunction
#' @importFrom checkmate assert check_choice check_class
setValidity("AlphaStableExtMO2FParam",
  function(object) {
    assert(combine = "and",
      check_choice(object@bf@scale, object@lambda),
      check_class(object@bf@original, "AlphaStableBernsteinFunction"))

    invisible(TRUE)
  })

#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions
#'   LinearBernsteinFunction PoissonBernsteinFunction
#' @importFrom checkmate assert check_choice check_class
setValidity("PoissonExtMO2FParam",
  function(object) {
    assert(combine = "and",
      check_choice(object@bf@scale, object@lambda),
      check_class(object@bf@original, "SumOfBernsteinFunctions"),
      check_class(object@bf@original@first, "LinearBernsteinFunction"),
      check_class(object@bf@original@second, "PoissonBernsteinFunction"))

      invisible(TRUE)
  })

#' @importFrom rmo ScaledBernsteinFunction SumOfBernsteinFunctions
#'   LinearBernsteinFunction ExponentialBernsteinFunction
#' @importFrom checkmate assert check_choice check_class
setValidity("ExponentialExtMO2FParam",
  function(object) {
    assert(combine = "and",
      check_choice(object@bf@scale, object@lambda),
      check_class(object@bf@original, "SumOfBernsteinFunctions"),
      check_class(object@bf@original@first, "LinearBernsteinFunction"),
      check_class(object@bf@original@second, "ExponentialBernsteinFunction"))

      invisible(TRUE)
  })

#' @importFrom checkmate qassert
setValidity("ExtGaussian2FParam",
  function(object) {
    qassert(object@lambda, "N1(0,)")
    qassert(object@nu, "N1[0,1]")

    invisible(TRUE)
  })

#' @importFrom copula getTheta
#' @importFrom checkmate qassert assert_true
setValidity("ExtArch2FParam",
  function(object) {
    qassert(object@lambda, "N1(0,)")
    qassert(object@nu, "N1")
    assert_true(object@dim == object@copula@dimension)
    assert_true(object@nu == getTheta(object@copula))

    invisible(TRUE)
  })
