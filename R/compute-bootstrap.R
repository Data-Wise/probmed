#' Parametric Bootstrap for P_med
#'
#' @keywords internal
.pmed_parametric_boot <- function(extract, x_ref, x_value, n_boot, ci_level, seed, ...) {
  if (!is.null(seed)) set.seed(seed)

  # Sample from multivariate normal
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' needed for parametric bootstrap. Install with: install.packages('MASS')")
  }

  theta_star <- MASS::mvrnorm(
    n = n_boot,
    mu = extract@estimates,
    Sigma = extract@vcov
  )

  # Extract indices for a, b, c'
  n_m <- length(extract@mediator_predictors)
  idx_a <- which(extract@mediator_predictors == extract@treatment)
  idx_b <- n_m + which(extract@outcome_predictors == extract@mediator)
  idx_c <- n_m + which(extract@outcome_predictors == extract@treatment)

  # Extract sigma values
  sigma_y <- extract@sigma_y
  sigma_m <- extract@sigma_m

  # Compute P_med for each bootstrap sample
  boot_estimates <- numeric(n_boot)

  for (i in seq_len(n_boot)) {
    boot_estimates[i] <- .pmed_core_simple(
      a = theta_star[i, idx_a],
      b = theta_star[i, idx_b],
      c_prime = theta_star[i, idx_c],
      x_ref = x_ref,
      x_value = x_value,
      sigma_y = sigma_y,
      sigma_m = sigma_m
    )
  }

  # Compute confidence interval
  alpha <- 1 - ci_level
  ci <- stats::quantile(boot_estimates, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)

  # Compute IE bootstrap estimates
  ie_boot_estimates <- theta_star[, idx_a] * theta_star[, idx_b]
  ie_ci <- stats::quantile(ie_boot_estimates, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)

  PmedResult(
    estimate = mean(boot_estimates, na.rm = TRUE),
    ci_lower = ci[1],
    ci_upper = ci[2],
    ci_level = ci_level,
    method = "parametric_bootstrap",
    n_boot = as.integer(n_boot),
    boot_estimates = boot_estimates,
    ie_estimate = mean(ie_boot_estimates, na.rm = TRUE),
    ie_ci_lower = ie_ci[1],
    ie_ci_upper = ie_ci[2],
    ie_boot_estimates = ie_boot_estimates,
    x_ref = x_ref,
    x_value = x_value,
    source_extract = extract,
    converged = TRUE
  )
}

#' Nonparametric Bootstrap for P_med
#'
#' @keywords internal
.pmed_nonparametric_boot <- function(extract, x_ref, x_value, n_boot, ci_level, seed, ...) {
  if (!is.null(seed)) set.seed(seed)

  data <- extract@data
  if (is.null(data)) {
    stop("Data is required for nonparametric bootstrap.")
  }
  n <- nrow(data)

  boot_estimates <- numeric(n_boot)
  ie_boot_estimates <- numeric(n_boot)

  for (i in seq_len(n_boot)) {
    # Resample data
    boot_indices <- sample(1:n, size = n, replace = TRUE)
    boot_data <- data[boot_indices, ]

    # Refit models
    # Refit models
    preds_m <- setdiff(extract@mediator_predictors, "(Intercept)")
    if (length(preds_m) == 0) preds_m <- "1"
    formula_m <- stats::reformulate(preds_m, response = extract@mediator)

    preds_y <- setdiff(extract@outcome_predictors, "(Intercept)")
    if (length(preds_y) == 0) preds_y <- "1"
    formula_y <- stats::reformulate(preds_y, response = extract@outcome)

    fit_m_boot <- stats::glm(formula_m, data = boot_data)
    fit_y_boot <- stats::glm(formula_y, data = boot_data)

    # Extract parameters
    a_boot <- stats::coef(fit_m_boot)[extract@treatment]
    b_boot <- stats::coef(fit_y_boot)[extract@mediator]
    c_boot <- stats::coef(fit_y_boot)[extract@treatment]
    sigma_m_boot <- stats::sigma(fit_m_boot)
    sigma_y_boot <- stats::sigma(fit_y_boot)

    # Compute P_med
    boot_estimates[i] <- .pmed_core_simple(
      a = a_boot,
      b = b_boot,
      c_prime = c_boot,
      x_ref = x_ref,
      x_value = x_value,
      sigma_y = sigma_y_boot,
      sigma_m = sigma_m_boot
    )

    # Compute IE
    ie_boot_estimates[i] <- a_boot * b_boot
  }

  # Compute confidence interval
  alpha <- 1 - ci_level
  ci <- stats::quantile(boot_estimates, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)

  # Compute IE CI
  ie_ci <- stats::quantile(ie_boot_estimates, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)

  PmedResult(
    estimate = mean(boot_estimates, na.rm = TRUE),
    ci_lower = ci[1],
    ci_upper = ci[2],
    ci_level = ci_level,
    method = "nonparametric_bootstrap",
    n_boot = as.integer(n_boot),
    boot_estimates = boot_estimates,
    ie_estimate = mean(ie_boot_estimates, na.rm = TRUE),
    ie_ci_lower = ie_ci[1],
    ie_ci_upper = ie_ci[2],
    ie_boot_estimates = ie_boot_estimates,
    x_ref = x_ref,
    x_value = x_value,
    source_extract = extract,
    converged = TRUE
  )
}
