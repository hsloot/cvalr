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
  function(.Object, ex_qmatrix) { # nolint
    if (!missing(ex_qmatrix)) {
      setExQMatrix(.Object) <- ex_qmatrix
      validObject(.Object)
    }

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
  definition = function(.Object, ex_intensities) { # nolint
    if (!missing(ex_intensities)) {
      setExIntensities(.Object) <- ex_intensities
      validObject(.Object)
    }

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
  definition = function(.Object, dim, bf) { # nolint
    if (!missing(dim) && !missing(bf)) {
      setDimension(.Object) <- dim
      setBernsteinFunction(.Object) <- bf
      validObject(.Object)
    }

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
setMethod("initialize", signature = "ExtMO2FParam", # nolint
  definition = function(.Object, # nolint
      dim, lambda, nu, rho = NULL, tau = NULL, alpha = NULL) {
    if (!missing(dim) && !missing(lambda) &&
          !(missing(nu) && missing(rho) && missing(tau) && missing(alpha))) {
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
    }

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
      dim, lambda, nu, rho = NULL, tau = NULL) {
    if (!missing(dim) && !missing(lambda) &&
          !(missing(nu) && missing(rho) && missing(tau))) {
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
    }

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
#' ExtArch2FParam(dim = 2L, lambda = 8e-2, rho = 4e-1, family = "Gumbel")
#' ExtArch2FParam(dim = 2L, lambda = 8e-2, tau = 4e-1, family = "Amh")
#' @importFrom copula archmCopula
setMethod("initialize", "ExtArch2FParam",
  function(.Object, # nolint
      dim = 2, lambda = 0.1, nu = 1, rho = NULL, tau = NULL,
      survival = TRUE,
      family = c("Clayton", "Frank", "Amh", "Gumbel", "Joe")) {
    family <- match.arg(family)
    .Object@survival <- survival
    .Object@copula <- archmCopula(family = family)
    if (!missing(dim) && !missing(lambda) &&
          !(missing(nu) && missing(rho) && missing(tau))) {
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
    }

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
#' ClaytonExtArch2FParam(dim = 2L, lambda = 8e-2, rho = 4e-1)
#' ClaytonExtArch2FParam(dim = 2L, lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "ClaytonExtArch2FParam",
  function(.Object, ..., survival = FALSE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Clayton"))
  })

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,FrankExtArch2FParam-method
#' @aliases initialize,FrankExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' FrankExtArch2FParam(dim = 2L, lambda = 8e-2, rho = 4e-1)
#' FrankExtArch2FParam(dim = 2L, lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "FrankExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Frank"))
  })

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,GumbelExtArch2FParam-method
#' @aliases initialize,GumbelExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' GumbelExtArch2FParam(dim = 2L, lambda = 8e-2, rho = 4e-1)
#' GumbelExtArch2FParam(dim = 2L, lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "GumbelExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Gumbel"))
  })

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,AmhExtArch2FParam-method
#' @aliases initialize,AmhExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' AmhExtArch2FParam(dim = 2L, lambda = 8e-2, rho = 4e-1)
#' AmhExtArch2FParam(dim = 2L, lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "AmhExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Amh"))
  })

#' @describeIn ExtArch2FParam-class Constructor
#' @aliases initialize,JoeExtArch2FParam-method
#' @aliases initialize,JoeExtArch2FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' JoeExtArch2FParam(dim = 2L, lambda = 8e-2, tau = 4e-1)
setMethod("initialize", "JoeExtArch2FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Joe"))
  })


#' @importFrom purrr map_lgl map_int reduce
#' @importFrom checkmate assert_true
setMethod("initialize", "H2ExMarkovParam",
  function(.Object, models = NULL, fraction = NULL) {
    if (!missing(models) && !missing(fraction)) {
      assert_true(all(map_lgl(models, is, class = "CalibrationParam")))
      dims <- map_int(models, getDimension)
      if (length(dims) > 1L) {
        partition <- reduce(
          dims[-1], ~{
            c(.x, list(max(c(0L, unlist(.x))) + 1:.y))
          }, .init = list())
      } else {
        partition <- list(1:dims)
      }

      .Object@dim <- dims[[1]]
      .Object@partition <- partition
      .Object@models <- models
      .Object@fraction <- fraction

      validObject(.Object)
    }

    invisible(.Object)
  })

#' @importFrom purrr imap
setMethod("initialize", "H2ExtMO3FParam",
  function(.Object,
    partition = list(1L:2L, 3L:5L), lambda = 1e-1, nu = c(0.2, 0.3),
    fraction = 0.5, rho = NULL, tau = NULL, alpha = NULL) {
  if (!missing(partition) && !missing(lambda) &&
        (!missing(nu) || !missing(rho) || !missing(tau) || !missing(alpha)) &&
        !missing(fraction)) {
    if (missing(nu)) {
      if (!is.null(rho)) {
        nu <- invRho(.Object, rho)
      } else if (!is.null(tau)) {
        nu <- invTau(.Object, tau)
      } else if (!is.null(alpha)) {
        nu <- invAlpha(.Object, alpha)
      }
    }

    dim <- length(unlist(partition))
    models <- imap(c(dim, map_int(partition, length)), ~{
      new(getModelName(.Object),
            dim = .x, lambda = lambda, nu = nu[[pmin(.y, 2L)]])
      })

    .Object@dim <- dim
    .Object@partition <- partition
    .Object@models <- models
    .Object@fraction <- fraction
    .Object@lambda <- lambda
    .Object@nu <- nu

    validObject(.Object)
  }

  invisible(.Object)
  })

#' @importFrom purrr imap
setMethod("initialize", "H2ExtGaussian3FParam",
  function(.Object,
      partition = list(1L:2L, 3L:5L), lambda = 1e-1, nu = c(0.2, 0.3),
      rho = NULL, tau = NULL) {
    if (!missing(partition) && !missing(lambda) &&
          (!missing(nu) || !missing(rho) || !missing(tau))) {
      if (missing(nu)) {
        if (!is.null(rho)) {
          nu <- invRho(.Object, rho)
        } else if (!is.null(tau)) {
          nu <- invTau(.Object, tau)
        }
      }

      dim <- length(unlist(partition))

      .Object@dim <- dim
      .Object@partition <- partition
      .Object@lambda <- lambda
      .Object@nu <- nu

      validObject(.Object)
    }

    invisible(.Object)
  })

#' @describeIn H2ExtArch3FParam-class Constructor
#' @aliases initialize,H2ExtArch3FParam-method
#' @aliases initialize,H2ExtArch3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param partition Partition of the components (only adjacent grouping allowed).
#' @param lambda Marginal intensity.
#' @param nu Dependency parameter (see [copula::archmCopula-class] and
#'   [copula::nacopula-class]).
#' @param rho Spearman's Rho.
#' @param tau Kendall's Tau.
#' @param survival Flag if survival copula is specified (default, except for Clayton)
#' @param family Name of the Archimedean copula family
#'   (see [copula::archmCopula-class]).
#'
#' @examples
#' H2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, tau = c(0.3, 0.5),
#'   survival = FALSE, family = "Clayton")
#' H2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, tau = c(0.3, 0.5),
#'   survival = TRUE, family = "Gumbel")
#'
#' @importFrom copula archmCopula onacopulaL iRho iTau
#' @importFrom purrr imap map
setMethod("initialize", "H2ExtArch3FParam",
  function(.Object,
      partition = list(1L:2L, 3L:5L),
      lambda = 1e-1, nu = c(0.2, 0.3), rho = NULL, tau = NULL,
      survival = TRUE,
      family = c("Clayton", "Frank", "AMH", "Gumbel", "Joe")) {
    family <- match.arg(family)
    .Object@family <- family
    .Object@survival <- survival
    acopula <- archmCopula(family = family)
    if (!missing(partition) && !missing(lambda) &&
          (!missing(nu) || !missing(rho) || !missing(tau))) {
      if (missing(nu)) {
        if (!is.null(rho)) {
          nu <- map_dbl(rho, ~copula::iRho(acopula, .))
        } else if (!is.null(tau)) {
          nu <- map_dbl(tau, ~copula::iTau(acopula, .))
        }
      }

      dim <- length(unlist(partition))

      .Object@dim <- as.integer(dim)
      .Object@partition <- partition
      .Object@lambda <- lambda
      .Object@nu <- nu

      .Object@copula <- onacopulaL(family = family, list(nu[[1]], c(), map(partition, ~list(nu[[2]], .))))

      validObject(.Object)
    }

    invisible(.Object)
  })

#' @describeIn H2ExtArch3FParam-class Constructor
#' @aliases initialize,ClaytonH2ExtArch3FParam-method
#' @aliases initialize,ClaytonH2ExtArch3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' ClaytonH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, tau = c(0.3, 0.5))
#' ClaytonH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, rho = c(0.3, 0.5))
setMethod("initialize", "ClaytonH2ExtArch3FParam",
  function(.Object, ..., survival = FALSE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Clayton"))
  })

#' @describeIn H2ExtArch3FParam-class Constructor
#' @aliases initialize,FrankH2ExtArch3FParam-method
#' @aliases initialize,FrankH2ExtArch3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' FrankH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, tau = c(0.3, 0.5))
#' FrankH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, rho = c(0.3, 0.5))
setMethod("initialize", "FrankH2ExtArch3FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Frank"))
  })

#' @describeIn H2ExtArch3FParam-class Constructor
#' @aliases initialize,GumbelH2ExtArch3FParam-method
#' @aliases initialize,GumbelH2ExtArch3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' GumbelH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, tau = c(0.3, 0.5))
#' GumbelH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, rho = c(0.3, 0.5))
setMethod("initialize", "GumbelH2ExtArch3FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Gumbel"))
  })

#' @describeIn H2ExtArch3FParam-class Constructor
#' @aliases initialize,AmhH2ExtArch3FParam-method
#' @aliases initialize,AmhH2ExtArch3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' AmhH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, tau = c(0.05, 0.2))
#' AmhH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, rho = c(0.05, 0.2))
setMethod("initialize", "AmhH2ExtArch3FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "AMH"))
  })

#' @describeIn H2ExtArch3FParam-class Constructor
#' @aliases initialize,JoeH2ExtArch3FParam-method
#' @aliases initialize,JoeH2ExtArch3FParam,ANY-method
#'
#' @inheritParams methods::initialize
#' @param ... Pass-through parameters.
#'
#' @examples
#' JoeH2ExtArch3FParam(
#'   partition = list(1:3, 4:6, 7:10, 11:15),
#'   lambda = 8e-2, tau = c(0.3, 0.5))
setMethod("initialize", "JoeH2ExtArch3FParam",
  function(.Object, ..., survival = TRUE) { # nolint
    invisible(callNextMethod(.Object, ..., survival = survival, family = "Joe"))
  })
