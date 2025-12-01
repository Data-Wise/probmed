#' S3 Class Wrappers for S7
#'
#' @description
#' S7 wrappers for S3 classes used in method dispatch
#'
#' @keywords internal
lm_class <- S7::new_S3_class("lm")

#' @keywords internal
glm_class <- S7::new_S3_class("glm")

#' Mediation Extract Base Class
#'
#' @description
#' S7 class for extracted mediation structures from fitted models.
#'
#' @param estimates Numeric vector of all model parameter estimates
#' @param mediator_predictors Character vector of mediator model predictor names
#' @param outcome_predictors Character vector of outcome model predictor names
#' @param a_path Numeric: treatment effect on mediator (a path)
#' @param b_path Numeric: mediator effect on outcome (b path)
#' @param c_prime Numeric: direct effect of treatment on outcome (c' path)
#' @param vcov Variance-covariance matrix of parameter estimates
#' @param sigma_m Numeric: residual standard deviation of mediator model
#' @param sigma_y Numeric: residual standard deviation of outcome model
#' @param data Data frame used for model fitting
#' @param n_obs Integer: number of observations
#' @param source_package Character: name of package that created the models
#' @param converged Logical: did the model(s) converge
#' @param treatment Character: name of treatment variable
#' @param mediator Character: name of mediator variable
#' @param outcome Character: name of outcome variable
#'
#' @export
MediationExtract <- S7::new_class(
  "MediationExtract",
  package = "probmed",
  properties = list(
    # Parameter estimates
    estimates = S7::class_numeric,

    # Model structure
    mediator_predictors = S7::class_character,
    outcome_predictors = S7::class_character,

    # Path coefficients
    a_path = S7::class_numeric,
    b_path = S7::class_numeric,
    c_prime = S7::class_numeric,

    # Variance-covariance matrix
    vcov = S7::class_any, # matrix class not exported by S7

    # Residual standard deviations
    sigma_m = S7::class_numeric,
    sigma_y = S7::class_numeric,

    # Data
    data = S7::new_property(
      class = S7::class_data.frame,
      default = data.frame()
    ),
    n_obs = S7::class_integer,

    # Metadata
    source_package = S7::class_character,
    converged = S7::class_logical,

    # Variable names
    treatment = S7::class_character,
    mediator = S7::class_character,
    outcome = S7::class_character
  ),
  validator = function(self) {
    if (length(self@a_path) != 1) {
      "a_path must be a single numeric value"
    } else if (length(self@b_path) != 1) {
      "b_path must be a single numeric value"
    } else if (self@n_obs < 1) {
      "n_obs must be positive"
    } else if (!self@converged) {
      warning("Model did not converge. Results may be unreliable.")
      NULL
    }
  }
)

#' P_med Result Class
#'
#' @description
#' S7 class for P_med computation results.
#'
#' @param estimate Numeric: P_med point estimate
#' @param ci_lower Numeric: lower bound of confidence interval
#' @param ci_upper Numeric: upper bound of confidence interval
#' @param ci_level Numeric: confidence level (e.g., 0.95)
#' @param method Character: inference method used
#' @param n_boot Integer: number of bootstrap samples (NA if not bootstrap)
#' @param boot_estimates Numeric vector: bootstrap distribution (empty if not bootstrap)
#' @param x_ref Numeric: reference treatment value
#' @param x_value Numeric: treatment value for contrast
#' @param source_extract MediationExtract object: source of the estimates
#' @param converged Logical: did computation converge
#' @param call Call object: original function call
#'
#' @export
PmedResult <- S7::new_class(
  "PmedResult",
  package = "probmed",
  properties = list(
    # Point estimate
    estimate = S7::class_numeric,

    # Confidence interval
    ci_lower = S7::class_numeric,
    ci_upper = S7::class_numeric,
    ci_level = S7::class_numeric,

    # Inference method
    method = S7::class_character,
    n_boot = S7::new_property(
      class = S7::class_integer,
      default = NA_integer_
    ),

    # Bootstrap distribution (if applicable)
    boot_estimates = S7::new_property(
      class = S7::class_numeric,
      default = numeric(0)
    ),

    # Treatment contrast
    x_ref = S7::class_numeric,
    x_value = S7::class_numeric,

    # Source information
    source_extract = S7::new_property(class = MediationExtract),

    # Computation details
    converged = S7::class_logical,
    call = S7::new_property(class = S7::class_any, default = NULL)
  ),
  validator = function(self) {
    if (self@estimate < 0 || self@estimate > 1) {
      warning("P_med outside [0,1]. Check model specification.")
    }

    if (!is.na(self@ci_lower) && !is.na(self@ci_upper)) {
      if (self@ci_lower > self@ci_upper) {
        "ci_lower must be <= ci_upper"
      }
    }

    if (!is.na(self@ci_level)) {
      if (self@ci_level <= 0 || self@ci_level >= 1) {
        "ci_level must be in (0, 1)"
      }
    }
  }
)
