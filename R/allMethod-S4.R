#' @include allClass-S4.R allGeneric-S4.R
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
