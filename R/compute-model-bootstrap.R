#' Model-based response bootstrap for P_med (opt-in)
#'
#' Holds the observed treatment and covariate rows fixed and, for each draw,
#' simulates a complete synthetic dataset from the fitted sequential models:
#' the mediator from its fitted conditional model, then the outcome from its
#' fitted model **using the simulated mediator**. Both equations are refitted
#' on those same rows (so coefficients and residual scales are re-estimated
#' jointly, keeping their pairing), and P_med is recomputed, standardized over
#' the original covariate rows when the model needs it.
#'
#' The point estimate is the fitted-model plug-in estimate (not the mean of
#' the draws, which is the legacy bootstraps' convention). The interval is the
#' percentile interval of the valid draws, reported only when at least 99% of
#' the requested draws are valid (finite estimate, converged supported fits,
#' no separation); otherwise both bounds are `NA` with status
#' `"bootstrap_incomplete"`. Failed draws are never recycled or replaced.
#'
#' Monte Carlo control: every draw evaluates P_med on one shared set of
#' uniforms `u` (common random numbers), so the simulation error is common to
#' all draws and does not widen the interval. That common error would instead
#' shift the whole interval, so each draw is re-centered as
#' `P_u(draw) - P_u(fit) + point`, where `P_u(fit)` is P_med at the fitted
#' coefficients on the same `u` (a control variate) and `point` uses
#' `100 * n_sim` simulated pairs. The remaining Monte Carlo error is the point
#' estimate's own, about `0.5 / sqrt(100 * n_sim)`.
#'
#' This models uncertainty under the specified sequential model and its
#' assumed confounding structure. It does not correct misspecification, and it
#' does not resample the covariate distribution.
#'
#' @keywords internal
.pmed_model_boot <- function(extract, x_ref, x_value, n_boot, ci_level, seed,
                             n_sim = 10000, ...) {
  if (is.null(extract@data)) {
    stop("method = \"model_bootstrap\" needs the source data.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  family_m <- if (is.null(extract@family_m)) stats::gaussian() else extract@family_m
  family_y <- if (is.null(extract@family_y)) stats::gaussian() else extract@family_y
  gaussian_m <- .is_gaussian(family_m)
  gaussian_y <- .is_gaussian(family_y)
  tx <- extract@treatment
  med <- extract@mediator

  design <- .design_matrices(extract)
  co <- .split_coefs(extract@estimates, design)
  n <- nrow(design$m)
  sigma_m <- extract@sigma_m
  sigma_y <- extract@sigma_y
  needs_std <- .needs_standardization(extract)

  # Point estimate: the fitted-model (standardized) plug-in.
  point <- .pmed_plugin(extract, x_ref, x_value, n_sim = 100L * n_sim)@estimate
  u <- .crn_uniforms(n_sim)
  p_u <- function(cm, cy, sm, sy) {
    std <- if (needs_std) {
      .std_linear_predictors(design, cm, cy, tx, med)
    } else list(m = NULL, y = NULL)
    .pmed_core_simple(
      a = cm[[tx]], b = cy[[med]],
      c_prime = if (tx %in% names(cy)) cy[[tx]] else 0,
      x_ref = x_ref, x_value = x_value,
      sigma_y = sy, sigma_m = sm, n_sim = n_sim,
      family_y = family_y, family_m = family_m,
      i_y = if ("(Intercept)" %in% names(cy)) cy[["(Intercept)"]] else 0,
      i_m = if ("(Intercept)" %in% names(cm)) cm[["(Intercept)"]] else 0,
      lp_m0 = std$m, lp_y0 = std$y, u = u
    )
  }
  # Control variate: the fitted model on the same u.
  anchor <- p_u(co$m, co$y,
                if (gaussian_m) sigma_m else NA_real_,
                if (gaussian_y) sigma_y else NA_real_)
  ie_point <- unname(co$m[[tx]] * co$y[[med]])

  eta_m <- drop(design$m %*% co$m)
  reasons <- character(n_boot)
  boot_estimates <- rep(NA_real_, n_boot)
  ie_boot <- rep(NA_real_, n_boot)

  for (i in seq_len(n_boot)) {
    # Simulate M | X, C, then Y | X, M*, C on the observed rows.
    m_star <- if (gaussian_m) {
      eta_m + stats::rnorm(n, 0, sigma_m)
    } else {
      stats::rbinom(n, 1, family_m$linkinv(eta_m))
    }
    dy <- design$y
    dy[, med] <- m_star
    eta_y <- drop(dy %*% co$y)
    y_star <- if (gaussian_y) {
      eta_y + stats::rnorm(n, 0, sigma_y)
    } else {
      stats::rbinom(n, 1, family_y$linkinv(eta_y))
    }

    if (!gaussian_m && .detect_separation(design$m, m_star) == "separated") {
      reasons[i] <- "separated"
      next
    }
    if (!gaussian_y && .detect_separation(dy, y_star) == "separated") {
      reasons[i] <- "separated"
      next
    }

    fit_m <- tryCatch(suppressWarnings(
      stats::glm.fit(design$m, m_star, family = family_m)), error = function(e) NULL)
    fit_y <- tryCatch(suppressWarnings(
      stats::glm.fit(dy, y_star, family = family_y)), error = function(e) NULL)
    if (is.null(fit_m) || is.null(fit_y)) {
      reasons[i] <- "fit_error"
      next
    }
    if (!isTRUE(fit_m$converged) || !isTRUE(fit_y$converged)) {
      reasons[i] <- "nonconvergence"
      next
    }
    cm <- fit_m$coefficients
    cy <- fit_y$coefficients
    if (!all(is.finite(cm)) || !all(is.finite(cy))) {
      reasons[i] <- "nonfinite"
      next
    }

    sm <- if (gaussian_m) sqrt(fit_m$deviance / fit_m$df.residual) else NA_real_
    sy <- if (gaussian_y) sqrt(fit_y$deviance / fit_y$df.residual) else NA_real_
    # Standardized over the ORIGINAL covariate rows (fixed target).
    est <- min(max(p_u(cm, cy, sm, sy) - anchor + point, 0), 1)
    if (!is.finite(est)) {
      reasons[i] <- "nonfinite"
      next
    }
    boot_estimates[i] <- est
    ie_boot[i] <- cm[[tx]] * cy[[med]]
    reasons[i] <- "ok"
  }

  valid <- reasons == "ok"
  n_valid <- sum(valid)
  complete <- n_valid >= ceiling(0.99 * n_boot)
  alpha <- 1 - ci_level
  probs <- c(alpha / 2, 1 - alpha / 2)
  if (complete) {
    ci <- stats::quantile(boot_estimates[valid], probs, names = FALSE)
    ie_ci <- stats::quantile(ie_boot[valid], probs, names = FALSE)
    status <- "ok"
    endpoint_status <- c(lower = "percentile", upper = "percentile")
  } else {
    ci <- ie_ci <- c(NA_real_, NA_real_)
    status <- "bootstrap_incomplete"
    endpoint_status <- c(lower = "failed", upper = "failed")
  }
  failures <- table(factor(reasons[!valid],
    levels = c("separated", "nonconvergence", "nonfinite", "fit_error")))

  PmedResult(
    estimate = point,
    ci_lower = ci[1],
    ci_upper = ci[2],
    ci_level = ci_level,
    method = "model_bootstrap",
    n_boot = as.integer(n_boot),
    boot_estimates = boot_estimates,
    ie_estimate = ie_point,
    ie_ci_lower = ie_ci[1],
    ie_ci_upper = ie_ci[2],
    ie_boot_estimates = ie_boot,
    x_ref = x_ref,
    x_value = x_value,
    source_extract = extract,
    converged = complete,
    diagnostics = list(
      interval_status = status,
      endpoint_status = endpoint_status,
      n_requested = as.integer(n_boot),
      n_valid = as.integer(n_valid),
      completion_rule = 0.99,
      failures = stats::setNames(as.integer(failures), names(failures)),
      n_sim = as.integer(n_sim)
    )
  )
}
