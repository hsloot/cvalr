#' @include allClass.R allGeneric.R
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


setMethod("invRho", "ExtMO2FParam",
  function(object, value) {
    invAlpha(object, 4 * value / (3 + value))
  })
setMethod("invRho", "ExtGaussian2FParam",
  function(object, value) {
    2 * sin(value * pi / 6)
  })
#' @importFrom copula iRho frankCopula
setMethod("invRho", "FrankExtArch2FParam",
  function(object, value) {
    copula::iRho(frankCopula(), value)
  })

setMethod("invTau", "ExtMO2FParam",
  function(object, value) {
    invAlpha(object, 2 * value / (1 + value))
  })
setMethod("invTau", "ExtGaussian2FParam",
  function(object, value) {
    sin(value * pi / 2)
  })
#' @importFrom copula iTau frankCopula
setMethod("invTau", "FrankExtArch2FParam",
  function(object, value) {
    copula::iTau(frankCopula(), value)
  })

setMethod("invAlpha", "CuadrasAugeExtMO2FParam",
  function(object, value) {
    value
  })
setMethod("invAlpha", "AlphaStableExtMO2FParam",
  function(object, value) {
    log2(2 - value)
  })
setMethod("invAlpha", "PoissonExtMO2FParam",
  function(object, value) {
    -log(1 - sqrt(value))
  })
setMethod("invAlpha", "ExponentialExtMO2FParam",
  function(object, value) {
    0.5 * (-3 + sqrt(1 + 8 / value))
  })
