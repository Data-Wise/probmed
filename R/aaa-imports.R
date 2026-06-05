#' @import S7
#' @importFrom stats coef vcov nobs glm predict rnorm quantile sd formula terms
#' @importFrom methods is
#' @keywords internal
"_PACKAGE"

# Package-level documentation
#' probmed: Probabilistic Effect Sizes for Mediation Analysis
#'
#' @description
#' Compute P_med, a scale-free probabilistic effect size for causal mediation
#' analysis. P_med = P(Y(x, M(x)) > Y(x, M(x*))) represents the probability that
#' the outcome under treatment with the mediator at its treated level exceeds
#' the outcome under treatment with the mediator at its control level (with a
#' one-half tie correction). The direct effect cancels, so P_med reflects
#' mediation only.
#'
#' @name probmed-package
#' @aliases probmed
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  # Register S7 methods for dynamic dispatch
  S7::methods_register()
  invisible()
}
