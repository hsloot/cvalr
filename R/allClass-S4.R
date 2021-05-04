#' @importFrom methods setClass setValidity setGeneric setMethod
#'   setReplaceMethod validObject is as new callNextMethod
NULL

#' Virtual super-class for calibration parameters
#'
#' A virtual superclass for all calibration parameters for multivariate
#' (portfolio) credit models.
#'
#' @slot dim The dimension (number of portfolio items)
#'
#' @details
#' It is assumed that the expected value of the credit derivative
#' payment streams can be represented as a function of *expected losses* and
#' other non-model dependent quantities.
#' *Expected losses* are the expected values of the average default counting
#' process \eqn{L} under a transformation \eqn{g} on times \eqn{t_1, \ldots, t_n},
#' i.e. \eqn{E[g(L_{t_1})], \ldots, E[g(L_{t_n})]}.
#'
#' A `CalibrationParam` object is supposed to be used to abstract the calculation
#' of these *expected losses* in functions which calculate the fair value of the
#' option with the asset pricing theorem.
#'
#' @seealso [ExMarkovParam-class] [ExMOParam-class] [ExtMOParam-class]
#'   [CuadrasAugeExtMO2FParam-class] [AlphaStableExtMO2FParam-class]
#'   [PoissonExtMO2FParam-class] [ExponentialExtMO2FParam-class]
#'
#' @export
setClass("CalibrationParam", # nolint
  contains = "VIRTUAL",
  slots = c(dim = "integer"))


#' Exchangeable Markovian calibration parameter
#'
#' Calibration parameter class for the general exchangeable model with a
#' Markovian *default counting process*.
#'
#' @slot qmatrix The \eqn{(d+1) \times (d+1)} Markov generator matrix of the
#'   default counting process.
#'
#' @details
#' The probability of \eqn{j > i} portfolio items being defaulted at time
#' \eqn{t > s} conditioned on \eqn{i} portfolio items being defaulted at time
#' \eqn{s} is
#'
#' \deqn{
#'   \mathbb{P}(Z_t = j \mid Z_s = i)
#'     = \delta_{i}^\top \operatorname{e}^{(t-s) Q} \delta_{j} .
#' }
#'
#' @export ExMarkovParam
ExMarkovParam <- setClass("ExMarkovParam", # nolint
  contains = "CalibrationParam",
  slots = c(ex_qmatrix = "matrix"))

#' Exchangeable Marshall--Olkin calibration parameter
#'
#' Calibration parameter classfor the general exchangeable model from the
#' Marshall--Olkin class.
#'
#' @slot ex_intensities The exchangeable intensities (see details)
#'
#' @details
#' The joint survival function of all portfolio items is assumed to be
#' \deqn{
#'   P(\tau > t)
#'     = \exp{(- a_0 t_{[1]} - \cdots - a_{d-1} t_{[d]})} ,
#' }
#' for \eqn{t_{[1]} \geq \cdots \geq t_{[d]}} begin the descendingly ordered
#' version of \eqn{t} and
#' \deqn{
#'   a_{i}
#'     = \sum_{l=0}^{d-i-1} \binom{d-i-1}{l} \lambda_{l+1} .
#' }
#'
#' @export ExMOParam
ExMOParam <- setClass("ExMOParam", # nolint
  contains = "ExMarkovParam",
  slots = c(ex_intensities = "numeric"))


#' Extendible Marshall--Olkin calibration parameter
#'
#' Calibration parameter class for the general extendible model from the
#' Marshall-Olkin class.
#'
#' @slot bf Bernstein function (see details)
#'
#' @details
#' The joint survival function of all portfolio items is assumed to be
#' \deqn{
#'   P(\tau > t)
#'     = \exp{(- a_{0} t_{[1]} - \cdots - a_{d-1} t_{[d]})} ,
#' }
#' for \eqn{t_{[1]} \geq \cdots \geq t_{[d]}} begin the descendingly ordered
#' version of \eqn{t} and
#' \deqn{
#'   a_{i}
#'     = \psi{(i+1)} - \psi{(i)} .
#' }
#'
#' @export ExtMOParam
ExtMOParam <- setClass("ExtMOParam", # nolint
  contains = "ExMOParam",
  slots = c(bf = "BernsteinFunction"))


#' Two-factor extendible Marshall--Olkin calibration parameter classes
#'
#' Calibration parameter classes with two parameters for the extendible
#' Marshall--Olkin family.
#'
#' @slot lambda The marginal rate
#' @slot nu Model specific dependence parameter
#'
#' @details
#' For all implemented families, the parameter `nu` can be replaced by
#' *Spearman's Rho* `rho`, *Kendall' Tau* `tau` or the
#' *(lower) tail dependence coefficient* `alpha`.
#' For all implemented families, the possible range for `rho`, `tau`, and
#' `alpha` is from zero to one (boundaries might not be included) and have a
#' one-to-one mapping to the model-specific parameter `nu`.
#' The link between lower tail dependence coefficient \eqn{\alpha} and
#' Spearman's Rho and Kendall's Tau is
#' \itemize{
#'   \item \eqn{\alpha = 4 \rho / (3 + \rho)} and \eqn{\rho = 3 \alpha / (4 - \alpha)}
#'   \item \eqn{\alpha = 2 \tau / (1 + \tau)} and \eqn{\tau = \alpha / (2 - \alpha)}
#' }
#'
#' @export
setClass("ExtMO2FParam", # nolint
  contains = c("ExtMOParam", "VIRTUAL"),
  slots = c(lambda = "numeric", nu = "numeric"))


#' @rdname ExtMO2FParam-class
#'
#' @section Cuadras-Augé calibration parameter class:
#' Corresponds to a Lévy subordinator which is a convex combination of
#' a pure-killing subordinator and a pure-drift subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = \nu + (1 - \nu) x}
#'   \item \eqn{\alpha = \nu}
#' }
#'
#' @export CuadrasAugeExtMO2FParam
CuadrasAugeExtMO2FParam <- setClass("CuadrasAugeExtMO2FParam", # nolint
  contains = "ExtMO2FParam")


#' @rdname ExtMO2FParam-class
#'
#' @section Alpha-stable calibration parameter class:
#' Corresponds to an \eqn{\alpha}-stable subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = x^\nu}
#'   \item \eqn{\nu = \log_2(2 - \alpha)} and \eqn{\alpha = 2 - 2^\nu}
#' }
#'
#' @export AlphaStableExtMO2FParam
AlphaStableExtMO2FParam <- setClass("AlphaStableExtMO2FParam", # nolint
  contains = "ExtMO2FParam")


#' @rdname ExtMO2FParam-class
#'
#' @section Poisson calibration parameter class:
#' Corresponds to a Lévy subrodinator which is a convex combination of a
#' Poisson subordinator with jump size `nu` and a pure-drift subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = \operatorname{e}^{-\nu}x + (1 - \operatorname{e}^{-x \nu})}
#'  \item \eqn{\nu = -log(1 - sqrt(\alpha))} and \eqn{\alpha = (1 - \operatorname{e}^{-\eta})}
#' }
#'
#' @export PoissonExtMO2FParam
PoissonExtMO2FParam <- setClass("PoissonExtMO2FParam", # nolint
  contains = "ExtMO2FParam")


#' @rdname ExtMO2FParam-class
#'
#' @section Exponential calibration parameter class:
#' Corresponds to a Lévy subordinator which is a convex combination of an
#' Exponential-jump compound Poisson process with rate `nu` and unit-intensity
#' and a pure-drift subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = (1 - 1 / (1 + \nu))x + 1 / (x + \nu)}
#'   \item \eqn{\nu = 0.5 \cdot (-3 + \sqrt{1 + 8 / \alpha})}
#'     and \eqn{\alpha = 2 / (1 + \nu) - 1 / (2 + \nu)}
#' }
#'
#' @export ExponentialExtMO2FParam
ExponentialExtMO2FParam <- setClass("ExponentialExtMO2FParam", # nolint
  contains = "ExtMO2FParam")


#' Two-factor extendible Gaussian calibration parameter classes
#'
#' Calibration parameter classes with two parameters for the extendible
#' Gaussian equi-correlation family with exponential margins.
#'
#' @slot lambda The marginal rate
#' @slot nu Model specific parameter (Pearson correlation)
#'
#' @details
#' For all implemented families, the parameter `nu` can be replaced by
#' *Spearman's Rho* `rho`, *Kendall' Tau* `tau`.
#' For all implemented families, the possible range for `rho` and `tau`
#' is from zero to one (boundaries might not be included) and have a
#' one-to-one mapping to the model-specific parameter `nu`.
#' The link between the Pearson correlation coefficient and
#' Spearman's Rho and Kendall's Tau is
#' \itemize{
#'   \item \eqn{\rho = 2 \sin(\rho_S \cdot \pi / 6)} and
#'     \eqn{\rho_S = 6 / \pi \cdot \arcsin(\rho/2)}
#'   \item \eqn{\rho = \sin(\tau \cdot \pi / 2)} and
#'     \eqn{\tau = 2 / \pi \cdot \arcsin(\rho)}
#' }
#'
#' @export ExtGaussian2FParam
ExtGaussian2FParam <- setClass("ExtGaussian2FParam", # nolint
  contains = "CalibrationParam",
  slots = c("lambda", "nu"))


#' Two-factor extendible Archimedean calibration parameter classes
#'
#' Calibration parameter classes with two parameters for the extendible
#' Archimedean families with exponential margins.
#'
#' @slot lambda The marginal rate
#' @slot nu Model specific parameter
#'
#' @details
#' For all implemented families, the parameter `nu` can be replaced by
#' *Spearman's Rho* `rho`, *Kendall's Tau* `tau`.
#' For all implemented families, the possible range for `rho` and `tau`
#' is from zero to one (boundaries might not be included) and have a
#' one-to-one mapping to the model-specific parameter `nu`.
#'
#' @importFrom copula iTau iRho tau rho frankCopula iPsi
#'
#' @export ExtArch2FParam
ExtArch2FParam <- setClass("ExtArch2FParam", # nolint
  contains = "CalibrationParam",
  slots = c(lambda = "numeric", nu = "numeric", survival = "logical",
    copula = "archmCopula"))


#' @rdname ExtArch2FParam-class
#'
#' @export ClaytonExtArch2FParam
ClaytonExtArch2FParam <- setClass("ClaytonExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "claytonCopula"))

#' @rdname ExtArch2FParam-class
#'
#' @export FrankExtArch2FParam
FrankExtArch2FParam <- setClass("FrankExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "frankCopula"))

#' @rdname ExtArch2FParam-class
#'
#' @export GumbelExtArch2FParam
GumbelExtArch2FParam <- setClass("GumbelExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "gumbelCopula"))

#' @rdname ExtArch2FParam-class
#'
#' @export AmhExtArch2FParam
AmhExtArch2FParam <- setClass("AmhExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "amhCopula"))

#' @rdname ExtArch2FParam-class
#'
#' @export JoeExtArch2FParam
JoeExtArch2FParam <- setClass("JoeExtArch2FParam", # nolint
  contains = "ExtArch2FParam",
  slots = c(lambda = "numeric", nu = "numeric", copula = "joeCopula"))


#' Virtual superclass for H2-exchangeable calibration parameters
#'
#' A virtual superclass for all 2-level hierarchically-exchangeable calibration
#' parameters for multivariate (portfolio) credit models.
#'
#' @slot partition Partition of the components (only adjacent grouping allowed)
#'
#' @export
setClass("H2ExCalibrationParam",
  contains = c("CalibrationParam", "VIRTUAL"),
  slots = c(partition = "list"))

#' H2-exchangeable Markovian calibration parameter
#'
#' Calibration parameter class for the general 2-level
#' hierarchically-exchangeable model with Markovian *default counting processes*
#' in the exchangeable sub-components.
#'
#' @slot models A list with the global and component models (of the type
#'   `ExMarkovParam-class`)
#' @slot fraction The proportion associated with the global model, see details.
#'
#' @details
#' We assume that \eqn{\tau} has the stochastic representation to be the
#' component-wise minimum of a global exchangeable Markovian-vector
#' \eqn{\tau^{(0)}} and a vector \eqn{(\tau^{(1)}, \ldots, \tau^{(J)})} with
#' independent exchangeable Markovian-vector sub-vectors \eqn{\tau^{(j)}}.
#'
#' @export H2ExMarkovParam
H2ExMarkovParam <- setClass("H2ExMarkovParam", # nolint
  contains = "H2ExCalibrationParam",
  slots = c(models = "list", fraction = "numeric"))

#' H2-Exchangeable Marshall--Olkin calibration parameter
#'
#' Calibration parameter class for the general 2-level
#' hierarchically-exchangeable model from the Marshall--Olkin class.
#'
#' @export H2ExMOParam
H2ExMOParam <- setClass("H2ExMOParam", # nolint
  contains = "H2ExMarkovParam")


#' H2-Extendible Marshall--Olkin calibration parameter
#'
#' Calibration parameter class for the general 2-level
#' hierarchically-etendiblew model from the Marshall--Olkin class.
#'
#' @export H2ExtMOParam
H2ExtMOParam <- setClass("H2ExtMOParam", # nolint
  contains = "H2ExMOParam")


#' Two-factor H2-extendible Marshall--Olkin calibration parameter classes
#'
#' Calibration parameter class with three parameters for the general 2-level
#' hierarchically-etendiblew model from the Marshall--Olkin class.
#'
#' @slot lambda The marginal rate
#' @slot nu Model-specific dependence parameters for the global- and the
#'   component-models.
#'
#' @details
#' For all implemented families, the parameters `nu`can be replaced by
#' *Spearman's Rho* `rho`, *Kendall's Tau* `tau`, or the *(lower) tail
#' dependence coefficient* `alpha`.
#' For all implemented families, the possible range for `rho`, `tau`, and
#' `alpha`, is between zero and one with the restriction that the global
#' parameter has to be smaller or equal than the corresponding component
#' parameter. Additionally, we support only that parameters of the same type are
#' provided, i.e. `rho`. The parameters have a one-to-one mapping to `nu`.
#'
#' @export
setClass("H2ExtMO3FParam", # nolint
  contains = c("H2ExtMOParam", "VIRTUAL"),
  slots = c(lambda = "numeric", nu = "numeric"))


#' @rdname H2ExtMO3FParam-class
#'
#' @section Cuadras-Augé calibration parameter class:
#' Corresponds to a Lévy subordinators which are a convex combination of
#' a pure-killing subordinator and a pure-drift subordinator.
#' \itemize{
#'   \item \eqn{\psi(x) = \nu + (1 - \nu) x}
#'   \item \eqn{\alpha = \nu}
#' }
#'
#' @export CuadrasAugeH2ExtMO3FParam
CuadrasAugeH2ExtMO3FParam <- setClass("CuadrasAugeH2ExtMO3FParam", # nolint
  contains = "H2ExtMO3FParam")


#' @rdname H2ExtMO3FParam-class
#'
#' @section Alpha-stable calibration parameter class:
#' Corresponds to \eqn{\alpha}-stable subordinators.
#' \itemize{
#'   \item \eqn{\psi(x) = x^\nu}
#'   \item \eqn{\nu = \log_2(2 - \alpha)} and \eqn{\alpha = 2 - 2^\nu}
#' }
#'
#' @export AlphaStableH2ExtMO3FParam
AlphaStableH2ExtMO3FParam <- setClass("AlphaStableH2ExtMO3FParam", # nolint
  contains = "H2ExtMO3FParam")


#' @rdname H2ExtMO3FParam-class
#'
#' @section Poisson calibration parameter class:
#' Corresponds to Lévy subordinators which are a convex combination of
#' Poisson subordinators with jump size `nu` and pure-drift subordinators.
#' \itemize{
#'   \item \eqn{\psi(x) = \operatorname{e}^{-\nu}x + (1 - \operatorname{e}^{-x \nu})}
#'  \item \eqn{\nu = -log(1 - sqrt(\alpha))} and \eqn{\alpha = (1 - \operatorname{e}^{-\eta})}
#' }
#'
#' @export PoissonH2ExtMO3FParam
PoissonH2ExtMO3FParam <- setClass("PoissonH2ExtMO3FParam", # nolint
  contains = "H2ExtMO3FParam")


#' @rdname H2ExtMO3FParam-class
#'
#' @section Exponential calibration parameter class:
#' Corresponds to Lévy subordinator which is a convex combination of
#' Exponential-jump compound Poisson processes with rate `nu` and unit-intensity
#' and a pure-drift subordinators.
#' \itemize{
#'   \item \eqn{\psi(x) = (1 - 1 / (1 + \nu))x + 1 / (x + \nu)}
#'   \item \eqn{\nu = 0.5 \cdot (-3 + \sqrt{1 + 8 / \alpha})}
#'     and \eqn{\alpha = 2 / (1 + \nu) - 1 / (2 + \nu)}
#' }
#'
#' @export ExponentialH2ExtMO3FParam
ExponentialH2ExtMO3FParam <- setClass("ExponentialH2ExtMO3FParam", # nolint
  contains = "H2ExtMO3FParam")


#' Three-factor H2-extendible Gaussian calibration parameter classes
#'
#' Calibration parameter classes with three parameters for the 2-level
#' hierarchically-extendible Gaussian family with exponential margins.
#'
#' @slot lambda The marginal rate
#' @slot nu Model-specific dependence parameters for the global- and the
#'   component-models.
#'
#' @details
#' For all implemented families, the parameters `nu` can be replaced by
#' *Spearman's Rho* `rho`, *Kendall's Tau* `tau`, or the *(lower) tail
#' dependence coefficient* `alpha`.
#' For all implemented families, the possible range for `rho`, `tau`, and
#' `alpha` is between zero and one with the restriction that the global
#' parameter has to be smaller or equal than the corresponding component
#' parameter. Additionally, we support only that parameters of the same type are
#' provided, i.e. `rho`. The parameters have a one-to-one mapping to `nu`.
#'
#' @export H2ExtGaussian3FParam
H2ExtGaussian3FParam <- setClass("H2ExtGaussian3FParam", # nolint
  contains = "H2ExCalibrationParam",
  slots = c(lambda = "numeric", nu = "numeric"))


#' Three-factor H2-extendible Archimedean calibration parameter classes
#'
#' Calibration parameter classes with three parameters for the 2-level
#' hierarchically-extendible Archimedean families with exponential margins.
#'
#' @slot lambda The marginal rate
#' @slot nu Model-specific dependence parameters for the global- and the
#'   component-models.
#' @slot partition Partition of the components (only adjacent grouping allowed)
#'
#' @details
#' For all implemented families, the parameters `nu` can be replaced by
#' *Spearman's Rho* `rho`, *Kendall's Tau* `tau`, or the *(lower) tail
#' dependence coefficient* `alpha`.
#' For all implemented families, the possible range for `rho`, `tau`, and
#' `alpha` is between zero and one with the restriction that the global
#' parameter has to be smaller or equal than the corresponding component
#' parameter. Additionally, we support only that parameters of the same type are
#' provided, i.e. `rho`. The parameters have a one-to-one mapping to `nu`.
#'
#' @export H2ExtArch3FParam
H2ExtArch3FParam <- setClass("H2ExtArch3FParam", # nolint
  contains = c("H2ExCalibrationParam"),
  slots = c(lambda = "numeric", nu = "numeric", family = "character",
    survival = "logical", copula = "outer_nacopula"))

#' @rdname H2ExtArch3FParam-class
#'
#' @export ClaytonH2ExtArch3FParam
ClaytonH2ExtArch3FParam <- setClass("ClaytonH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")

#' @rdname H2ExtArch3FParam-class
#'
#' @export FrankH2ExtArch3FParam
FrankH2ExtArch3FParam <- setClass("FrankH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")

#' @rdname H2ExtArch3FParam-class
#'
#' @export GumbelH2ExtArch3FParam
GumbelH2ExtArch3FParam <- setClass("GumbelH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")

#' @rdname H2ExtArch3FParam-class
#'
#' @export AmhH2ExtArch3FParam
AmhH2ExtArch3FParam <- setClass("AmhH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")

#' @rdname H2ExtArch3FParam-class
#'
#' @export JoeH2ExtArch3FParam
JoeH2ExtArch3FParam <- setClass("JoeH2ExtArch3FParam", # nolint
  contains = "H2ExtArch3FParam")
