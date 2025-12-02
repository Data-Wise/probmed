#' Extract Mediation Structure
#'
#' @description
#' Generic function to extract mediation structure (a, b, c' paths and
#' variance-covariance matrices) from fitted models. Supports multiple
#' model types including lm/glm, lavaan SEM objects, and mediation package
#' objects.
#'
#' @param object Fitted model object (lm, glm, lavaan, or mediate)
#' @param ... Additional arguments passed to methods. Common arguments include:
#'   \itemize{
#'     \item \code{treatment}: Character name of treatment variable
#'     \item \code{mediator}: Character name of mediator variable
#'     \item \code{outcome}: Character name of outcome variable (auto-detected for lavaan)
#'     \item \code{model_y}: Outcome model (for lm/glm methods)
#'     \item \code{data}: Data frame (for lm/glm methods)
#'   }
#'
#' @return \code{MediationExtract} object containing:
#'   \itemize{
#'     \item Path coefficients (a, b, c')
#'     \item Variance-covariance matrix
#'     \item Residual standard deviations
#'     \item Variable names and metadata
#'   }
#'
#' @examples
#' # Example 1: Extract from lm objects
#' set.seed(123)
#' data <- data.frame(X = rnorm(100), C = rnorm(100))
#' data$M <- 0.5 * data$X + 0.3 * data$C + rnorm(100)
#' data$Y <- 0.4 * data$M + 0.2 * data$X + rnorm(100)
#'
#' model_m <- lm(M ~ X + C, data = data)
#' model_y <- lm(Y ~ M + X + C, data = data)
#' extract <- extract_mediation(model_m, model_y = model_y,
#'                               treatment = "X", mediator = "M", data = data)
#' print(extract)
#'
#' \donttest{
#' # Example 2: Extract from lavaan
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#'   library(lavaan)
#'   model <- 'M ~ a*X + C
#'             Y ~ b*M + cp*X + C'
#'   fit <- sem(model, data = data)
#'   extract <- extract_mediation(fit, treatment = "X", mediator = "M")
#'   print(extract)
#' }
#'
#' # Example 3: Extract from mediate object
#' if (requireNamespace("mediation", quietly = TRUE)) {
#'   library(mediation)
#'   med_out <- mediation::mediate(model_m, model_y,
#'                                  treat = "X", mediator = "M", boot = FALSE)
#'   extract <- extract_mediation(med_out)
#'   print(extract)
#' }
#' }
#'
#' @seealso \code{\link{pmed}} for computing P_med from extracted objects
#' @export
extract_mediation <- S7::new_generic(
  "extract_mediation",
  dispatch_args = "object"
)

#' Compute P_med: Probabilistic Effect Size for Mediation
#'
#' @description
#' Compute \eqn{P_{med}}, a scale-free probabilistic effect size for mediation
#' analysis, along with the traditional Indirect Effect (\eqn{a \times b}).
#' Provides point estimates and bootstrap confidence intervals.
#'
#' \eqn{P_{med}} represents \eqn{P(Y_{x^*, M_x} > Y_{x, M_x})}, the probability
#' that the counterfactual outcome under control with mediator at treated level
#' exceeds the outcome under treatment.
#'
#' @param object Either:
#'   \itemize{
#'     \item A \code{formula} for the outcome model (most common)
#'     \item A \code{MediationExtract} object from \code{extract_mediation()}
#'   }
#' @param ... Additional arguments passed to methods (see method documentation)
#'
#' @return \code{PmedResult} object containing:
#'   \itemize{
#'     \item \code{estimate}: P_med point estimate
#'     \item \code{ci_lower}, \code{ci_upper}: Confidence interval bounds
#'     \item \code{ie_estimate}: Indirect Effect point estimate
#'     \item \code{ie_ci_lower}, \code{ie_ci_upper}: IE confidence interval
#'     \item \code{boot_estimates}: Bootstrap distribution (if applicable)
#'     \item \code{method}: Inference method used
#'   }
#'
#' @section Methods:
#' Available inference methods via \code{method} argument:
#' \itemize{
#'   \item \code{"parametric_bootstrap"} (default): Fast, assumes normality
#'   \item \code{"nonparametric_bootstrap"}: Robust to assumptions, slower
#'   \item \code{"plugin"}: Point estimate only, no CI
#' }
#'
#' @examples
#' # Basic example with formula interface
#' set.seed(123)
#' n <- 200
#' data <- data.frame(X = rnorm(n), C = rnorm(n))
#' data$M <- 0.5 * data$X + 0.3 * data$C + rnorm(n)
#' data$Y <- 0.4 * data$M + 0.2 * data$X + 0.2 * data$C + rnorm(n)
#'
#' # Compute with parametric bootstrap
#' result <- pmed(
#'   Y ~ X + M + C,
#'   formula_m = M ~ X + C,
#'   data = data,
#'   treatment = "X",
#'   mediator = "M",
#'   method = "parametric_bootstrap",
#'   n_boot = 500
#' )
#' print(result)
#'
#' \dontrun{
#' # View bootstrap distribution
#' summary(result)
#' plot(result)
#' }
#'
#' @seealso \code{\link{extract_mediation}} for extracting from fitted models
#' @export
pmed <- S7::new_generic("pmed", dispatch_args = "object")
