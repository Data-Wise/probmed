#' Extract from GLM/LM
#'
#' @name extract_mediation-lm
S7::method(extract_mediation, lm_class) <- function(object,
                                                    treatment,
                                                    mediator,
                                                    model_y = NULL,
                                                    data = NULL,
                                                    ...) {
  # Validate inputs
  if (is.null(model_y)) {
    stop(
      "For lm/glm objects, provide both models:\n",
      "  model_m: mediator model (M ~ X + C)\n",
      "  model_y: outcome model (Y ~ X + M + C)\n",
      "Use: extract_mediation(model_m, model_y = model_y, ...)"
    )
  }

  # Extract from models (object is model_m)
  model_m <- object

  # Extract coefficients early for validation
  coefs_m <- stats::coef(model_m)
  coefs_y <- stats::coef(model_y)

  # Validate treatment variable exists in mediator model
  if (!(treatment %in% names(coefs_m))) {
    stop(
      "Treatment variable '", treatment, "' not found in mediator model.\n",
      "Available predictors: ", paste(names(coefs_m), collapse = ", ")
    )
  }

  # Validate mediator variable exists in outcome model
  if (!(mediator %in% names(coefs_y))) {
    stop(
      "Mediator variable '", mediator, "' not found in outcome model.\n",
      "Available predictors: ", paste(names(coefs_y), collapse = ", ")
    )
  }

  # Validate treatment variable exists in outcome model
  if (!(treatment %in% names(coefs_y))) {
    stop(
      "Treatment variable '", treatment, "' not found in outcome model.\n",
      "For mediation analysis, the outcome model should include both treatment and mediator.\n",
      "Available predictors: ", paste(names(coefs_y), collapse = ", ")
    )
  }

  # Get data
  if (is.null(data)) {
    data <- model_m$model
    if (is.null(data)) {
      warning(
        "No data available in model object. ",
        "Bootstrap methods will not work. ",
        "Provide data explicitly via the 'data' argument."
      )
    }
  }

  # Get path coefficients
  a_path <- coefs_m[treatment]
  b_path <- coefs_y[mediator]
  c_prime <- coefs_y[treatment]

  # Combined parameter vector
  estimates <- c(coefs_m, coefs_y)

  # Combined vcov (block diagonal)
  vcov_m <- stats::vcov(model_m)
  vcov_y <- stats::vcov(model_y)
  n_m <- length(coefs_m)
  n_y <- length(coefs_y)
  vcov_combined <- matrix(0, n_m + n_y, n_m + n_y)
  vcov_combined[1:n_m, 1:n_m] <- vcov_m
  vcov_combined[(n_m + 1):(n_m + n_y), (n_m + 1):(n_m + n_y)] <- vcov_y

  # Extract residual standard deviations (only for Gaussian models)
  # For GLMs with non-Gaussian families (binomial, poisson, etc.),
  # residual variance is not meaningful in the same way

  # Check mediator model family
  family_m <- if (inherits(model_m, "glm")) {
    model_m$family$family
  } else {
    "gaussian" # lm objects are Gaussian
  }

  # Check outcome model family
  family_y <- if (inherits(model_y, "glm")) {
    model_y$family$family
  } else {
    "gaussian"
  }

  # Extract sigma only for Gaussian families
  sigma_m <- if (family_m == "gaussian") {
    stats::sigma(model_m)
  } else {
    NA_real_ # Not applicable for non-Gaussian
  }

  sigma_y <- if (family_y == "gaussian") {
    stats::sigma(model_y)
  } else {
    NA_real_
  }

  # Create MediationExtract
  MediationExtract(
    estimates = estimates,
    mediator_predictors = names(coefs_m),
    outcome_predictors = names(coefs_y),
    a_path = unname(a_path),
    b_path = unname(b_path),
    c_prime = unname(c_prime),
    vcov = vcov_combined,
    sigma_m = sigma_m,
    sigma_y = sigma_y,
    data = data,
    n_obs = stats::nobs(model_m),
    source_package = "stats",
    converged = model_m$converged %||% TRUE,
    treatment = treatment,
    mediator = mediator,
    outcome = all.vars(stats::formula(model_y))[1]
  )
}
