#' Mediation Extract Base Class
#'
#' @description
#' S7 class for extracted mediation structures from fitted models.
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
    vcov = S7::class_matrix,
    
    # Data
    data = S7::class_data.frame,
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
    
    if (self@ci_level <= 0 || self@ci_level >= 1) {
      "ci_level must be in (0, 1)"
    }
  }
)