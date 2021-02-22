#' @include allClass-S4.R allGeneric-S4.R
NULL

setMethod("initialize", "ExMarkovParam",
  function(.Object, # nolint
      qmatrix = 0.1 * 0.5 * matrix(c(-3, 0, 0, 2, -2, 0, 1, 2, 0), nrow=3, ncol=3)) {
    setQMatrix(.Object) <- qmatrix
    validObject(.Object)

    invisible(.Object)
  })

setMethod("initialize", "ExMOParam",
  definition = function(.Object, # nolint
                        ex_intensities = 0.1 * 0.5 * c(1, 1)) {
    setExIntensities(.Object) <- ex_intensities
    validObject(.Object)

    invisible(.Object)
  })

setMethod("initialize", "ExtMOParam",
  definition = function(.Object, # nolint
      dim = 2L,
      bf = ScaledBernsteinFunction(
        scale = 0.05,
        original = SumOfBernsteinFunctions(
          first = ConstantBernsteinFunction(constant = 0.5),
          second = LinearBernsteinFunction(scale = 0.5)
        )
      )) {
    setDimension(.Object) <- dim
    setBernsteinFunction(.Object) <- bf
    validObject(.Object)

    invisible(.Object)
  })

setMethod("initialize", signature = "ExtMO2FParam",
  definition = function(.Object, # nolint
      dim = 2, lambda = 0.1, nu = 0.5, rho = NULL, tau = NULL, alpha = NULL) {
    if (missing(nu)) {
      if (!is.null(rho)) {
        nu <- invRho(.Object, rho)
      } else if (!is.null(tau)) {
        nu <- invTau(.Object, tau)
      } else if (!is.null(alpha)) {
        nu <- invAlpha(.Object, alpha)
      }
    }

    setDimension(.Object) <- dim
    setBernsteinFunction(.Object) <- constructBernsteinFunction(.Object, lambda, nu)
    validObject(.Object)

    invisible(.Object)
  })

setMethod("initialize", signature = "ExtGaussian2FParam",
  definition = function(.Object, # nolint
      dim = 2, lambda = 0.1, nu = 0.5, rho = NULL, tau = NULL) {
    if (missing(nu)) {
      if (!is.null(rho)) {
        nu <- invRho(.Object, rho)
      } else if (!is.null(tau)) {
        nu <- invTau(.Object, tau)
      }
    }

    setDimension(.Object) <- dim
    setLambda(.Object) <- lambda
    setNu(.Object) <- nu
    validObject(.Object)

    invisible(.Object)
  })

setMethod("initialize", "FrankExtArch2FParam",
  function(.Object, # nolint
      dim = 2, lambda = 0.1, nu = 0.5, rho = NULL, tau = NULL) {
    if (missing(nu)) {
      if (!is.null(rho)) {
        nu <- invRho(.Object, rho)
      } else if (!is.null(tau)) {
        nu <- invTau(.Object, tau)
      }
    }

    setDimension(.Object) <- dim
    setLambda(.Object) <- lambda
    setNu(.Object) <- nu
    validObject(.Object)

    invisible(.Object)
  })
