#' @include allClass.R allGeneric.R
NULL

#' Probability vector and expected value for calibration parameter
#'
#' Calculates the probability vector and expected values of the average
#' default counting process \eqn{L}.
#'
#' @param object The calibration parameter object
#' @param t Point-in-time
#'
#' @docType methods
#' @export
setGeneric("probability_vector",
  function(object, t) {
    standardGeneric("probability_vector")
  })

#' @rdname probability_vector
#' @aliases probability_vector,ExMarkovParam
#'
#' @examples
#' probability_vector(CuadrasAugeExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#' probability_vector(AlphaStableExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#' probability_vector(PoissonExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#' probability_vector(ExponentialExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#'
#' @importFrom expm expm
#' @export
setMethod("probability_vector", "ExMarkovParam",
  function(object, t) {
    as.vector(expm(t * object@qmatrix)[1, ])
  })

#' @rdname probability_vector
#' @aliases probability_vector,ExtGaussian2FParam
#'
#' @examples
#' probability_vector(ExtGaussian2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#'
#' @importFrom stats integrate pexp pnorm dnorm qnorm
#' @export
setMethod("probability_vector", "ExtGaussian2FParam",
  function(object, t) {
    sapply(0:object@dim,
      function(k) {
        int_res <- integrate(
          function(x) {
            ldp <- pnorm(
              (qnorm(pexp(t, rate = object@lambda)) - sqrt(object@nu) * x) /
              (sqrt(1 - object@nu)),
              log.p=TRUE, lower.tail = TRUE
            )
            lsp <- pnorm(
              (qnorm(pexp(t, rate = object@lambda)) - sqrt(object@nu) * x) /
              (sqrt(1 - object@nu)),
              log.p=TRUE, lower.tail = FALSE
            )
            sapply(
              exp(k * ldp + (object@dim-k) * lsp) * dnorm(x),
              function(v) {
                multiply_binomial_coefficient(v, object@dim, k)
              })
          }, lower = -Inf, upper = Inf, rel.tol = .Machine$double.eps^0.5
        )

        int_res$value
      })
  })

#' @rdname probability_vector
#' @aliases probability_vector,FrankExtArch2FParam
#'
#' @examples
#' probability_vector(FrankExtArch2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#'
#' @importFrom stats pexp
#' @importFrom copula iPsi frankCopula
#' @export
setMethod("probability_vector", "FrankExtArch2FParam",
  function(object, t) {
    tt <- copula::iPsi(frankCopula(object@nu),
                        pexp(t, rate = object@lambda))
    dn <- function(m) {
      pexp(object@nu)^m / (object@nu * m)
    }
    n <- max(which(dn(2 ^ (1:10)) > .Machine$double.eps))
    sapply(0:object@dim,
      function(k) {
        multiply_binomial_coefficient(
          sum(
            pexp(tt, rate = (1:n), lower.tail = FALSE) ^ k *
              pexp(tt, rate = (1:n), lower.tail = TRUE) ^ (object@dim - k) *
              dn((1:n))
          ),
          object@dim, k)
      })
  })
