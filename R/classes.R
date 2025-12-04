#' S3 Class Wrappers for S7
#'
#' @description
#' S7 wrappers for S3 classes used in method dispatch
#'
#' @keywords internal
lm_class <- S7::new_S3_class("lm")

#' @keywords internal
glm_class <- S7::new_S3_class("glm")

#' @keywords internal
lavaan_class <- S7::new_S3_class("lavaan") # Placeholder if needed, but we use dynamic S4

#' @keywords internal
mediate_class <- S7::new_S3_class("mediate")


# NOTE: MediationExtract has been replaced by medfit::MediationData
# probmed now depends on medfit for the core mediation data structure.
# This allows sharing infrastructure across the mediationverse ecosystem.

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
#' @param ie_estimate Numeric: Indirect Effect (NIE) point estimate
#' @param ie_ci_lower Numeric: lower bound of NIE confidence interval
#' @param ie_ci_upper Numeric: upper bound of NIE confidence interval
#' @param ie_boot_estimates Numeric vector: bootstrap distribution of NIE
#' @param x_ref Numeric: reference treatment value
#' @param x_value Numeric: treatment value for contrast
#' @param source_extract medfit::MediationData object: source of the estimates
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

    # Indirect Effect (NIE)
    ie_estimate = S7::new_property(class = S7::class_numeric, default = NA_real_),
    ie_ci_lower = S7::new_property(class = S7::class_numeric, default = NA_real_),
    ie_ci_upper = S7::new_property(class = S7::class_numeric, default = NA_real_),
    ie_boot_estimates = S7::new_property(
      class = S7::class_numeric,
      default = numeric(0)
    ),

    # Treatment contrast
    x_ref = S7::class_numeric,
    x_value = S7::class_numeric,

    # Source information (medfit::MediationData object)
    source_extract = S7::class_any,

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
