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
    qassert(object@survival, "B1")
    assert_true(object@dim == object@copula@dimension)
    assert_true(object@nu == getTheta(object@copula))

    invisible(TRUE)
  })


#' @importFrom purrr map_lgl
#' @importFrom checkmate qassert qtest assert_true
setValidity("H2ExCalibrationParam",
  function(object) {
    assert_true(all(map_lgl(object@partition, ~{
      qtest(.x, "I+(0,)")
      })))
    assert_true(all(1L:object@dim == unlist(object@partition)))

    invisible(TRUE)
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

#' @importFrom checkmate qassert
setValidity("H2ExtMO3FParam",
  function(object) {
    qassert(object@lambda, "N1(0,)")
    qassert(object@nu, "N2(0,)")

    invisible(TRUE)
  })

#' @importFrom checkmate qassert
setValidity("H2ExtGaussian3FParam",
  function(object) {
    qassert(object@lambda, "N1(0,)")
    qassert(object@nu, "N2(0,)")

    invisible(TRUE)
  })

#' @importFrom checkmate qassert
setValidity("H2ExtArch3FParam",
  function(object) {
    qassert(object@lambda, "N1(0,)")
    qassert(object@nu, "N2(0,)")
    qassert(object@survival, "B1")
    assert_true(object@dim == dim(object@copula))
    assert_true(check_equal(object@nu, c(object@copula@copula@theta, object@copula@childCops[[1]]@copula@theta)))

    invisible(TRUE)
  })
