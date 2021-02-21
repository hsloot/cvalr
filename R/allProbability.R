#' @include allClass.R allGeneric.R
NULL

#' @describeIn CalibrationParam-class
#'   returns the probability vector for the average default count process \eqn{L}.
#'
#' @param object The calibration parameter object
#' @param times Point-in-time
#'
#' @export
setGeneric("probability_distribution",
  function(object, times) {
    standardGeneric("probability_distribution")
  })

#' @describeIn ExMarkovParam-class
#'   returns the probability vector for the average default count process \eqn{L}.
#' @aliases probability_distribution,ExMarkovParam-method
#'
#' @inheritParams probability_distribution
#'
#' @examples
#' probability_distribution(CuadrasAugeExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#' probability_distribution(AlphaStableExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#' probability_distribution(PoissonExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#' probability_distribution(ExponentialExtMO2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#'
#' @importFrom expm expm
#' @importFrom checkmate qassert
#' @export
setMethod("probability_distribution", "ExMarkovParam",
  function(object, times) {
    qassert(times, "N+[0,)")
    sapply(times, function(t) expm(t * object@qmatrix)[1, ])
  })

#' @describeIn ExtGaussian2FParam-class
#'   returns the probability vector for the average default count process \eqn{L}.
#' @aliases probability_distribution,ExtGaussian2FParam-method
#'
#' @inheritParams probability_distribution
#'
#' @examples
#' probability_distribution(ExtGaussian2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#'
#' @importFrom stats integrate pexp pnorm dnorm qnorm
#' @importFrom checkmate qassert
#' @export
setMethod("probability_distribution", "ExtGaussian2FParam",
  function(object, times) {
    qassert(times, "N+[0,)")
    times <- qnorm(pexp(times, rate = object@lambda))
    out <- simplify2array(outer(0:object@dim, times,
      Vectorize(function(k, t) {
        int_res <- integrate(
          function(x) {
            ldp <- pnorm(
              (t - sqrt(object@nu) * x) / (sqrt(1 - object@nu)),
              log.p=TRUE, lower.tail = TRUE
            )
            lsp <- pnorm(
              (t - sqrt(object@nu) * x) / (sqrt(1 - object@nu)),
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
      })), higher=FALSE)

    if (1L == nrow(out) || 1L == ncol(out))
      out <- as.vector(out)

    out
  })

#' @describeIn ExtArch2FParam-class
#'   returns the probability vector for the average default count process \eqn{L}.
#' @aliases probability_distribution,FrankExtArch2FParam-method
#'
#' @inheritParams probability_distribution
#'
#' @examples
#' probability_distribution(FrankExtArch2FParam(
#'   dim = 50, lambda = 0.05, rho = 0.4), 0.3)
#'
#' @importFrom stats pexp
#' @importFrom copula iPsi frankCopula
#' @importFrom checkmate qassert
#' @export
setMethod("probability_distribution", "FrankExtArch2FParam",
  function(object, times) {
    qassert(times, "N+[0,)")
    times <- copula::iPsi(frankCopula(object@nu),
                        pexp(times, rate = object@lambda))
    dn <- function(m) {
      pexp(object@nu)^m / (object@nu * m)
    }
    n <- max(which(dn(2 ^ (1:10)) > .Machine$double.eps))
    out <- simplify2array(outer(0:object@dim, times,
      Vectorize(function(k, t) {
        multiply_binomial_coefficient(
          sum(
            pexp(t, rate = (1:n), lower.tail = FALSE) ^ k *
              pexp(t, rate = (1:n), lower.tail = TRUE) ^ (object@dim - k) *
              dn(1:n)
          ),
          object@dim, k)
      })), higher=FALSE)

    if (1L == nrow(out) || 1L == ncol(out))
      out <- as.vector(out)

    out
  })
