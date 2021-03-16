#' @include allClass-S4.R allGeneric-S4.R
NULL

#' @describeIn ExMarkovParam-class Constructor
#' @aliases initialize,ExMarkovParam-method
#' @aliases initialize,ExMarkovParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ex_qmatrix (Exchangeable) Q-matrix of the default counting process.
#'
#' @examples
#' ExMarkovParam(
#'  ex_qmatrix = matrix(
#'    c(-0.07647059, 0, 0, 0.05294118, -0.05, 0, 0.02352941, 0.05, 0),
#'    nrow = 3L, ncol = 3L))
setMethod("initialize", "ExMarkovParam",
  function(.Object, # nolint
      ex_qmatrix = 0.1 * 0.5 * matrix(c(-3, 0, 0, 2, -2, 0, 1, 2, 0), nrow=3, ncol=3)) {
    setExQMatrix(.Object) <- ex_qmatrix
    validObject(.Object)

    invisible(.Object)
  })

#' @describeIn ExMOParam-class Constructor
#' @aliases initialize,ExMOParam-method
#' @aliases initialize,ExMOParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ex_intensities (Scaled) exchangeable intensities of the exchangeable
#'   Marshall-Olkin distribution.
#'
#' @examples
#' ExMOParam(ex_intensities = c(0.02647059, 0.02352941))
setMethod("initialize", "ExMOParam",
  definition = function(.Object, # nolint
                        ex_intensities = 0.1 * 0.5 * c(1, 1)) {
    setExIntensities(.Object) <- ex_intensities
    validObject(.Object)

    invisible(.Object)
  })

#' @describeIn ExtMOParam-class Constructor
#' @aliases initialize,ExtMOParam-method
#' @aliases initialize,ExtMOParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param dim Dimension.
#' @param bf Bernstein function.
#'
#' @examples
#' ExtMOParam(
#'   dim = 2,
#'   bf = rmo::ScaledBernsteinFunction(
#'     scale = 0.05,
#'     original = rmo::SumOfBernsteinFunctions(
#'       first = rmo::ConstantBernsteinFunction(constant = 0.4),
#'       second = rmo::LinearBernsteinFunction(scale = 1 - 0.4))
#'     ))
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

#' @describeIn ExtMO2FParam-class Constructor
#' @aliases initialize,ExtMO2FParam-method
#' @aliases initialize,ExtMO2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param dim Dimension.
#' @param lambda Marginal intensity.
#' @param nu Dependence parameter.
#' @param rho Spearman's Rho.
#' @param tau Kendall's Tau.
#' @param alpha Bivariate lower tail dependence coefficient
#'
#' @examples
#' CuadrasAugeExtMO2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
#' AlphaStableExtMO2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
#' PoissonExtMO2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
#' ExponentialExtMO2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
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

#' @describeIn ExtGaussian2FParam-class Constructor
#' @aliases initialize,ExtGaussian2FParam-method
#' @aliases initialize,ExtGaussian2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param dim Dimension.
#' @param lambda Marginal intensity.
#' @param nu Dependence parameter.
#' @param rho Spearman's Rho.
#' @param tau Kendall's Tau.
#'
#' @examples
#' ExtGaussian2FParam(dim = 2L, lambda = 0.05, rho = 0.4)
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

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,ExtArch2FParam-method
#' @aliases initialize,ExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param dim Dimension.
#' @param lambda Marginal intensity.
#' @param nu Dependence parameter (see [copula::archmCopula-class]).
#' @param rho Spearman's Rho.
#' @param tau Kendall's Tau.
#' @param survival Flag if survival copula is specified (default, except for Clayton)
#' @param family Name of the Archimedean copula family
#'   (see [copula::archmCopula-class]).
#'
#' @examples
#' ExtArch2FParam()
#' ExtArch2FParam(lambda = 8e-2, rho = 4e-1, family = "gumbel")
#' ExtArch2FParam(lambda = 8e-2, tau = 4e-1, family = "amh")
#' @importFrom copula archmCopula
setMethod("initialize", "ExtArch2FParam",
  function(.Object, # nolint
      dim = 2, lambda = 0.1, nu = 1, rho = NULL, tau = NULL,
      survival = TRUE,
      family = c("clayton", "frank", "amh", "gumbel", "joe")) {
    family <- match.arg(family)
    .Object@survival <- survival
    .Object@copula <- archmCopula(family = family)
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

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,ClaytonExtArch2FParam-method
#' @aliases initialize,ClaytonExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' ClaytonExtArch2FParam()
#' ClaytonExtArch2FParam(lambda = 8e-2, rho = 4e-1)
#' ClaytonExtArch2FParam(lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "ClaytonExtArch2FParam",
  function(.Object, ..., survival = FALSE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "clayton"))
  })

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,FrankExtArch2FParam-method
#' @aliases initialize,FrankExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' FrankExtArch2FParam()
#' FrankExtArch2FParam(lambda = 8e-2, rho = 4e-1)
#' FrankExtArch2FParam(lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "FrankExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "frank"))
  })

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,GumbelExtArch2FParam-method
#' @aliases initialize,GumbelExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' GumbelExtArch2FParam()
#' GumbelExtArch2FParam(lambda = 8e-2, rho = 4e-1)
#' GumbelExtArch2FParam(lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "GumbelExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "gumbel"))
  })

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,AmhExtArch2FParam-method
#' @aliases initialize,AmhExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' AmhExtArch2FParam()
#' AmhExtArch2FParam(lambda = 8e-2, rho = 4e-1)
#' AmhExtArch2FParam(lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "AmhExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "amh"))
  })

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,JoeExtArch2FParam-method
#' @aliases initialize,JoeExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' JoeExtArch2FParam()
#' JoeExtArch2FParam(lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "JoeExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "joe"))
  })
