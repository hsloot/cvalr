#' @include allClass-S4.R 
NULL

#' @importFrom rmo SumOfBernsteinFunctions LinearBernsteinFunction
#'   ConstantBernsteinFunction
setMethod("constructBernsteinFunction", "CuadrasAugeExtMO2FParam",
  function(object, lambda, nu, ...) {
    ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = 1 - nu),
        second = ConstantBernsteinFunction(constant = nu))
    )
  })

#' @importFrom rmo SumOfBernsteinFunctions AlphaStableBernsteinFunction
setMethod("constructBernsteinFunction", "AlphaStableExtMO2FParam",
  function(object, lambda, nu, ...) {
    ScaledBernsteinFunction(
      scale = lambda,
      original = AlphaStableBernsteinFunction(alpha = nu)
    )
  })

#' @importFrom rmo SumOfBernsteinFunctions LinearBernsteinFunction
#'   PoissonBernsteinFunction
setMethod("constructBernsteinFunction", "PoissonExtMO2FParam",
  function(object, lambda, nu, ...) {
    ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = exp(-nu)),
        second = PoissonBernsteinFunction(lambda = 1, eta = nu)
      )
    )
  })

#' @importFrom rmo SumOfBernsteinFunctions LinearBernsteinFunction
#'   ExponentialBernsteinFunction
setMethod("constructBernsteinFunction", "ExponentialExtMO2FParam",
  function(object, lambda, nu, ...) {
    ScaledBernsteinFunction(
      scale = lambda,
      original = SumOfBernsteinFunctions(
        first = LinearBernsteinFunction(scale = 1 - 1 / (1 + nu)),
        second = ExponentialBernsteinFunction(lambda = nu)
      )
    )
  })


#' @importFrom checkmate qassert
setMethod("invRho", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    invAlpha(object, 4 * value / (3 + value))
  })
#' @importFrom checkmate qassert
setMethod("invRho", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    2 * sin(value * pi / 6)
  })
#' @importFrom copula iRho frankCopula
#' @importFrom checkmate qassert
setMethod("invRho", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    copula::iRho(object@copula, value)
  })

#' @importFrom checkmate qassert
setMethod("invTau", "ExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    invAlpha(object, 2 * value / (1 + value))
  })
#' @importFrom checkmate qassert
setMethod("invTau", "ExtGaussian2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    sin(value * pi / 2)
  })
#' @importFrom copula iTau frankCopula
#' @importFrom checkmate qassert
setMethod("invTau", "ExtArch2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    copula::iTau(object@copula, value)
  })

#' @importFrom checkmate qassert
setMethod("invAlpha", "CuadrasAugeExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    value
  })
#' @importFrom checkmate qassert
setMethod("invAlpha", "AlphaStableExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    log2(2 - value)
  })
#' @importFrom checkmate qassert
setMethod("invAlpha", "PoissonExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    -log(1 - sqrt(value))
  })
#' @importFrom checkmate qassert
setMethod("invAlpha", "ExponentialExtMO2FParam",
  function(object, value) {
    qassert(value, "N1[0,1]")
    0.5 * (-3 + sqrt(1 + 8 / value))
  })

setMethod("getModelName", "H2ExMarkovParam",
  function(object) {
    "ExMarkovParam"
  })
setMethod("getModelName", "H2ExMOParam",
  function(object) {
    "ExMOParam"
  })
setMethod("getModelName", "H2ExtMOParam",
  function(object) {
    "ExtMOParam"
  })
setMethod("getModelName", "H2ExtMO3FParam",
  function(object) {
    "ExtMO2FParam"
  })
setMethod("getModelName", "CuadrasAugeH2ExtMO3FParam",
  function(object) {
    "CuadrasAugeExtMO2FParam"
  })
setMethod("getModelName", "AlphaStableH2ExtMO3FParam",
  function(object) {
    "AlphaStableExtMO2FParam"
  })
setMethod("getModelName", "PoissonH2ExtMO3FParam",
  function(object) {
    "PoissonExtMO2FParam"
  })
setMethod("getModelName", "ExponentialH2ExtMO3FParam",
  function(object) {
    "ExponentialExtMO2FParam"
  })
setMethod("getModelName", "H2ExtGaussian3FParam",
  function(object) {
    "ExtGaussian2FParam"
  })
setMethod("getModelName", "H2ExtArch3FParam",
  function(object) {
    "ExtArch2FParam"
  })
setMethod("getModelName", "ClaytonH2ExtArch3FParam",
  function(object) {
    "ClaytonExtArch2FParam"
  })
setMethod("getModelName", "FrankH2ExtArch3FParam",
  function(object) {
    "FrankExtArch2FParam"
  })
setMethod("getModelName", "GumbelH2ExtArch3FParam",
  function(object) {
    "GumbelExtArch2FParam"
  })
setMethod("getModelName", "AmhH2ExtArch3FParam",
  function(object) {
    "AmhExtArch2FParam"
  })
setMethod("getModelName", "JoeH2ExtArch3FParam",
  function(object) {
    "JoeExtArch2FParam"
  })

#' @importFrom checkmate qassert
setMethod("invRho", "H2ExtMO3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    invAlpha(object, 4 * value / (3 + value))
  })
#' @importFrom checkmate qassert
setMethod("invRho", "H2ExtGaussian3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    2 * sin(value * pi / 6)
  })
#' @importFrom copula iRho
#' @importFrom purrr map_dbl
#' @importFrom checkmate qassert
setMethod("invRho", "H2ExtArch3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    map_dbl(value, ~copula::iRho(object@copula@copula, .))
  })

#' @importFrom checkmate qassert
setMethod("invTau", "H2ExtMO3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    invAlpha(object, 2 * value / (1 + value))
  })
#' @importFrom checkmate qassert
setMethod("invTau", "H2ExtGaussian3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    sin(value * pi / 2)
  })
#' @importFrom copula iTau
#' @importFrom checkmate qassert
setMethod("invTau", "H2ExtArch3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    copula::iTau(object@copula@copula, value)
  })

#' @include utils.R
#' @importFrom checkmate qassert
setMethod("invAlpha", "CuadrasAugeH2ExtMO3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    adjacent_differences(value) / c(object@fraction, 1 - object@fraction)
  })
#' @include utils.R
#' @importFrom checkmate qassert
setMethod("invAlpha", "AlphaStableH2ExtMO3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    value <- adjacent_differences(value) / c(object@fraction, 1 - object@fraction)
    log2(2 - value)

  })
#' @include utils.R
#' @importFrom checkmate qassert
setMethod("invAlpha", "PoissonH2ExtMO3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    value <- adjacent_differences(value) / c(object@fraction, 1 - object@fraction)
    -log(1 - sqrt(value))
  })
#' @include utils.R
#' @importFrom checkmate qassert
setMethod("invAlpha", "ExponentialH2ExtMO3FParam",
  function(object, value) {
    qassert(value, "N2[0,1]")
    value <- adjacent_differences(value) / c(object@fraction, 1 - object@fraction)
    0.5 * (-3 + sqrt(1 + 8 / value))
  })
