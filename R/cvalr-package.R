## usethis namespace: start
#' @useDynLib cvalr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
#' @importFrom magrittr %>%
NULL

if(getRversion() >= "2.15.1") utils::globalVariables(c(".")) # nolint

#' cvalr: Credit derivative valuation in R
#'
#' The package contains a framework to valuate credit derivatives (at the moment
#' portfolio CDS swaps and CDO tranches) for various multivariate models (e.g.
#' exchangeable and hierarchical Marshall--Olkin / Archimedean / Gaussian
#' models).
#'
#' @section Calibration parameter:
#' A * calibration parameter* represents a family of multivariate models.  or
#' each model the method `expected_value()` can be used to calculate arbitrary
#' expectations on the *default counting process*. The implementation depends on
#' the model and can involve Monte-Carlo simulation. For specific derivatives,
#' the methods `expected_*_equation`, e.g. [expected_pcds_equation()] and
#' [expected_cdo_equation()], should be used to allow for optimized calculations of
#' the expectation (in particular for the Gaussian case).
#' - At the moment, all models have exponential margins. This is not a
#'   limitation for derivatives which are not path dependent (payoffs only
#'   depend on state values at the payoff date), as the time can be transformed
#'   to account for other margins.
#' - Most exchangeable models are implemented with two factors `lambda` and
#'   `nu`, where the former is the marginal rate and the latter is a dependence
#'   parameter which is able to represent a range between zero and one for all
#'   classical bivariate dependence measures. The methods `setAlpha()`,
#'   `setRho()`, and `setTau()` can be used to set the dependence parameter such
#'   that the model has a certain *(lower) tail-correlation coefficient*,
#'   *Spearman's Rho*, or *Kendall's Tau*.
#'
#' @docType package
#' @name cvalr
NULL
