#' Compute P_med from Formula
#'
#' @name pmed-formula
#' @param object Formula for outcome model (Y ~ X + M + C)
#' @param formula_m Formula for mediator model (M ~ X + C)
#' @param data Data frame
#' @param treatment Character: treatment variable name
#' @param mediator Character: mediator variable name
#' @param family_y Family for outcome model (default: gaussian())
#' @param family_m Family for mediator model (default: gaussian())
#' @param x_ref Reference treatment value (default: 0)
#' @param x_value Treatment value (default: 1)
#' @param method Inference method: "parametric_bootstrap", "nonparametric_bootstrap", "plugin"
#' @param n_boot Number of bootstrap samples (default: 1000)
#' @param ci_level Confidence level (default: 0.95)
#' @param seed Random seed for reproducibility
#' @param ... Additional arguments
#'
#' @return PmedResult object
#'
#' @examples
#' # Toy example: Simple mediation model
#' # Generate data where X affects Y through M
#' set.seed(123)
#' n <- 100
#' data <- data.frame(
#'   X = rnorm(n),
#'   C = rnorm(n)
#' )
#' data$M <- 0.5 * data$X + 0.3 * data$C + rnorm(n, sd = 0.5)
#' data$Y <- 0.4 * data$M + 0.2 * data$X + 0.2 * data$C + rnorm(n, sd = 0.5)
#'
#' # Compute P_med using plugin estimator (fast, no CI)
#' result_plugin <- pmed(
#'   Y ~ X + M + C,
#'   formula_m = M ~ X + C,
#'   data = data,
#'   treatment = "X",
#'   mediator = "M",
#'   method = "plugin"
#' )
#' print(result_plugin)
#'
#' \donttest{
#' # With parametric bootstrap for confidence intervals
#' result_boot <- pmed(
#'   Y ~ X + M + C,
#'   formula_m = M ~ X + C,
#'   data = data,
#'   treatment = "X",
#'   mediator = "M",
#'   method = "parametric_bootstrap",
#'   n_boot = 200,  # Use more (e.g., 1000+) in practice
#'   seed = 456
#' )
#' print(result_boot)
#' }
S7::method(pmed, S7::class_formula) <- function(object,  # formula_y
                                                formula_m,
                                                data,
                                                treatment,
                                                mediator,
                                                family_y = stats::gaussian(),
                                                family_m = stats::gaussian(),
                                                x_ref = 0,
                                                x_value = 1,
                                                method = c("parametric_bootstrap",
                                                          "nonparametric_bootstrap",
                                                          "plugin"),
                                                n_boot = 1000,
                                                ci_level = 0.95,
                                                seed = NULL,
                                                ...) {
  
  method <- match.arg(method)
  
  # Store call
  call <- match.call()
  
  # Fit models
  fit_m <- stats::glm(formula_m, data = data, family = family_m)
  fit_y <- stats::glm(object, data = data, family = family_y)
  
  # Extract structure
  extract <- extract_mediation(
    fit_m,
    treatment = treatment,
    mediator = mediator,
    model_y = fit_y,
    data = data
  )
  
  # Compute P_med
  result <- .pmed_compute(
    extract = extract,
    x_ref = x_ref,
    x_value = x_value,
    method = method,
    n_boot = n_boot,
    ci_level = ci_level,
    seed = seed,
    ...
  )
  
  # Add call
  result@call <- call
  
  return(result)
}

#' Compute P_med from MediationExtract
#'
#' @name pmed-MediationExtract
S7::method(pmed, MediationExtract) <- function(object,
                                                x_ref = 0,
                                                x_value = 1,
                                                method = c("parametric_bootstrap",
                                                          "nonparametric_bootstrap",
                                                          "plugin"),
                                                n_boot = 1000,
                                                ci_level = 0.95,
                                                seed = NULL,
                                                ...) {
  
  method <- match.arg(method)
  
  .pmed_compute(
    extract = object,
    x_ref = x_ref,
    x_value = x_value,
    method = method,
    n_boot = n_boot,
    ci_level = ci_level,
    seed = seed,
    ...
  )
}