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

  # Resolve coefficient columns BY NAME. medfit names the joint estimate vector
  # m_<pred> / y_<pred> (see extract_mediation), so these are unambiguous and
  # robust to column ordering (positional indexing previously selected the
  # intercepts by mistake).
  nm_a <- paste0("m_", extract@treatment)
  nm_b <- paste0("y_", extract@mediator)
  nm_c <- paste0("y_", extract@treatment)
  a_s <- theta_star[, nm_a]
  b_s <- theta_star[, nm_b]
  c_s <- theta_star[, nm_c]
  cn <- colnames(theta_star)
  iy_s <- if ("y_(Intercept)" %in% cn) theta_star[, "y_(Intercept)"] else rep(0, n_boot)
  im_s <- if ("m_(Intercept)" %in% cn) theta_star[, "m_(Intercept)"] else rep(0, n_boot)

  # Extract sigma + family (held fixed across parametric replicates)
  sigma_y <- extract@sigma_y
  sigma_m <- extract@sigma_m
  family_y <- extract@family_y
  family_m <- extract@family_m

  # Compute P_med for each bootstrap sample
  boot_estimates <- numeric(n_boot)

  for (i in seq_len(n_boot)) {
    boot_estimates[i] <- .pmed_core_simple(
      a = a_s[i],
      b = b_s[i],
      c_prime = c_s[i],
      x_ref = x_ref,
      x_value = x_value,
      sigma_y = sigma_y,
      sigma_m = sigma_m,
      family_y = family_y,
      family_m = family_m,
      i_y = iy_s[i],
      i_m = im_s[i]
    )
  }

  # Compute confidence interval
  alpha <- 1 - ci_level
  ci <- stats::quantile(boot_estimates, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)

  # Compute IE bootstrap estimates
  ie_boot_estimates <- a_s * b_s
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

  # Preserve the original model families so refits use the correct link
  # (default gaussian when absent).
  family_m <- if (is.null(extract@family_m)) stats::gaussian() else extract@family_m
  family_y <- if (is.null(extract@family_y)) stats::gaussian() else extract@family_y

  # Build model formulas once
  preds_m <- setdiff(extract@mediator_predictors, "(Intercept)")
  if (length(preds_m) == 0) preds_m <- "1"
  formula_m <- stats::reformulate(preds_m, response = extract@mediator)

  preds_y <- setdiff(extract@outcome_predictors, "(Intercept)")
  if (length(preds_y) == 0) preds_y <- "1"
  formula_y <- stats::reformulate(preds_y, response = extract@outcome)

  boot_estimates <- numeric(n_boot)
  ie_boot_estimates <- numeric(n_boot)

  for (i in seq_len(n_boot)) {
    # Resample data
    boot_indices <- sample(1:n, size = n, replace = TRUE)
    boot_data <- data[boot_indices, ]

    # Refit models on the correct family/link
    fit_m_boot <- stats::glm(formula_m, data = boot_data, family = family_m)
    fit_y_boot <- stats::glm(formula_y, data = boot_data, family = family_y)

    cm <- stats::coef(fit_m_boot)
    cy <- stats::coef(fit_y_boot)

    # Extract parameters (residual SD only meaningful for Gaussian)
    a_boot <- cm[extract@treatment]
    b_boot <- cy[extract@mediator]
    c_boot <- cy[extract@treatment]
    im_boot <- if ("(Intercept)" %in% names(cm)) cm[["(Intercept)"]] else 0
    iy_boot <- if ("(Intercept)" %in% names(cy)) cy[["(Intercept)"]] else 0
    sigma_m_boot <- if (.is_gaussian(family_m)) stats::sigma(fit_m_boot) else NA_real_
    sigma_y_boot <- if (.is_gaussian(family_y)) stats::sigma(fit_y_boot) else NA_real_

    # Compute P_med
    boot_estimates[i] <- .pmed_core_simple(
      a = a_boot,
      b = b_boot,
      c_prime = c_boot,
      x_ref = x_ref,
      x_value = x_value,
      sigma_y = sigma_y_boot,
      sigma_m = sigma_m_boot,
      family_y = family_y,
      family_m = family_m,
      i_y = iy_boot,
      i_m = im_boot
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
