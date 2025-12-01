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
#' analysis. P_med = P(Y_{x*,M_x} > Y_{x,M_x}) represents the probability that
#' the outcome under treatment with the mediator at its control level exceeds
#' the outcome under control.
#'
#' @docType package
#' @name probmed-package
NULL

.onLoad <- function(libname, pkgname) {
  # Register S7 classes
  invisible()
}