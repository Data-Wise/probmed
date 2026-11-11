#' Core P_med Computation Dispatcher
#'
#' @keywords internal
.pmed_compute <- function(extract, x_ref, x_value, method, n_boot, ci_level, seed,
                          n_sim = 10000, ...) {
  # Reject unsupported families/structures before any random draw.
  .check_supported_extract(extract, method = method)
  # A separated binary fit has no MLE; stop rather than simulate from it.
  separation <- .check_original_separation(extract)
  result <- switch(method,
    parametric_bootstrap = .pmed_parametric_boot(
      extract, x_ref, x_value, n_boot, ci_level, seed, n_sim = n_sim, ...
    ),
    nonparametric_bootstrap = .pmed_nonparametric_boot(
      extract, x_ref, x_value, n_boot, ci_level, seed, n_sim = n_sim, ...
    ),
    plugin = .pmed_plugin(
      extract, x_ref, x_value, n_sim = n_sim, seed = seed, ...
    ),
    mbco = .pmed_mbco(
      extract, x_ref, x_value, ci_level = ci_level, ...
    ),
    model_bootstrap = .pmed_model_boot(
      extract, x_ref, x_value, n_boot, ci_level, seed, n_sim = n_sim, ...
    )
  )
  if (length(separation)) {
    diagnostics <- result@diagnostics
    diagnostics$separation <- separation
    result@diagnostics <- diagnostics
  }
  result
}

#' Plugin Estimator (Point Estimate Only)
#'
#' @keywords internal
.pmed_plugin <- function(extract, x_ref, x_value, n_sim = 10000, seed = NULL, ...) {
  # Match bootstrap semantics: explicit seeds reset and then advance caller RNG.
  if (!is.null(seed)) set.seed(seed)
  # Extract parameters
  a <- extract@a_path
  b <- extract@b_path
  c_prime <- extract@c_prime
  sigma_y <- extract@sigma_y
  sigma_m <- extract@sigma_m
  std <- .std_for(extract, extract@estimates)

  # Compute P_med
  pmed_est <- .pmed_core_simple(
    a = a,
    b = b,
    c_prime = c_prime,
    x_ref = x_ref,
    x_value = x_value,
    sigma_y = sigma_y,
    sigma_m = sigma_m,
    n_sim = n_sim,
    family_y = extract@family_y,
    family_m = extract@family_m,
    i_y = .named_or(extract@estimates, "y_(Intercept)"),
    i_m = .named_or(extract@estimates, "m_(Intercept)"),
    lp_m0 = std$m,
    lp_y0 = std$y
  )

  # Compute Indirect Effect (a * b)
  ie_est <- a * b

  PmedResult(
    estimate = pmed_est,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    ci_level = NA_real_,
    method = "plugin",
    n_boot = NA_integer_,
    ie_estimate = ie_est,
    ie_ci_lower = NA_real_,
    ie_ci_upper = NA_real_,
    x_ref = x_ref,
    x_value = x_value,
    source_extract = extract,
    converged = TRUE
  )
}

#' Look up a named coefficient, returning a default when absent
#'
#' @keywords internal
.named_or <- function(estimates, name, default = 0) {
  if (!is.null(estimates) && name %in% names(estimates)) {
    unname(estimates[[name]])
  } else {
    default
  }
}

#' Is a family Gaussian (NULL is treated as Gaussian)?
#'
#' @keywords internal
.is_gaussian <- function(family) {
  is.null(family) || identical(family$family, "gaussian")
}

#' Simple P_med Core Computation
#'
#' Computes the mediation estimand
#' \eqn{P_{med} = P(Y(x, M(x)) > Y(x, M(x^*))) + \tfrac12 P(Y(x, M(x)) = Y(x, M(x^*)))}:
#' the mediator is drawn under BOTH treatment levels independently while the
#' treatment is held at `x_value` for both outcomes, so the direct effect
#' `c_prime` cancels and the result reflects mediation only. The intercepts and
#' families matter only for non-Gaussian models (they cancel in the Gaussian
#' linear contrast); their defaults keep Gaussian behaviour unchanged.
#' `lp_m0` / `lp_y0`, when given, are per-row covariate linear predictors
#' (see `.std_linear_predictors()`); each pair then draws one analysis row and
#' uses it in both arms, standardizing over the observed covariates.
#' `u`, when given, is an `n_sim` x 5 matrix of uniforms (see `.crn_uniforms()`)
#' used in place of fresh RNG draws via the inverse CDF. Reusing one `u` across
#' bootstrap draws (common random numbers) makes their Monte Carlo error move
#' together, so it no longer inflates the bootstrap spread.
#'
#' @keywords internal
.pmed_core_simple <- function(a, b, c_prime, x_ref, x_value,
                              sigma_y = 1, sigma_m = 1, n_sim = 10000,
                              family_y = NULL, family_m = NULL,
                              i_y = 0, i_m = 0,
                              lp_m0 = NULL, lp_y0 = NULL, u = NULL) {
  # Residual SDs are NULL/NA for non-Gaussian models; treat as 0 (unused there,
  # since the family branch draws responses directly).
  if (is.null(sigma_m) || is.na(sigma_m)) sigma_m <- 0
  if (is.null(sigma_y) || is.na(sigma_y)) sigma_y <- 0

  gaussian_m <- .is_gaussian(family_m)
  gaussian_y <- .is_gaussian(family_y)

  # --- Draw the mediator under BOTH treatment levels, independently ---
  # E[M | X] on the link scale; intercept cancels for Gaussian but is needed
  # for a nonlinear link.
  if (!is.null(lp_m0)) {
    # Standardize over the analysis rows: one row per pair, shared by both arms.
    rows <- if (is.null(u)) {
      sample.int(length(lp_m0), n_sim, replace = TRUE)
    } else {
      pmin(floor(u[, 1] * length(lp_m0)) + 1L, length(lp_m0))
    }
    i_m <- lp_m0[rows]
    i_y <- lp_y0[rows]
  }
  mu_m_t <- i_m + a * x_value # mediator mean under treatment
  mu_m_c <- i_m + a * x_ref # mediator mean under control
  if (!is.null(u)) {
    if (gaussian_m) {
      m_t <- mu_m_t + sigma_m * stats::qnorm(u[, 2])
      m_c <- mu_m_c + sigma_m * stats::qnorm(u[, 3])
    } else {
      m_t <- as.integer(u[, 2] < family_m$linkinv(mu_m_t))
      m_c <- as.integer(u[, 3] < family_m$linkinv(mu_m_c))
    }
  } else if (gaussian_m) {
    m_t <- stats::rnorm(n_sim, mean = mu_m_t, sd = sigma_m)
    m_c <- stats::rnorm(n_sim, mean = mu_m_c, sd = sigma_m)
  } else {
    # Non-Gaussian mediator (e.g. binary): draw responses on the correct scale
    m_t <- stats::rbinom(n_sim, 1, family_m$linkinv(mu_m_t))
    m_c <- stats::rbinom(n_sim, 1, family_m$linkinv(mu_m_c))
  }

  # --- Outcomes with treatment held at x_value for BOTH (c' cancels) ---
  eta_t <- i_y + b * m_t + c_prime * x_value # outcome lin. predictor, treated mediator
  eta_c <- i_y + b * m_c + c_prime * x_value # outcome lin. predictor, control mediator
  if (!is.null(u)) {
    if (gaussian_y) {
      y_t <- eta_t + sigma_y * stats::qnorm(u[, 4])
      y_c <- eta_c + sigma_y * stats::qnorm(u[, 5])
    } else {
      y_t <- as.integer(u[, 4] < family_y$linkinv(eta_t))
      y_c <- as.integer(u[, 5] < family_y$linkinv(eta_c))
    }
  } else if (gaussian_y) {
    y_t <- eta_t + stats::rnorm(n_sim, 0, sigma_y)
    y_c <- eta_c + stats::rnorm(n_sim, 0, sigma_y)
  } else {
    # Non-Gaussian outcome (e.g. binary): map through the link, draw responses
    y_t <- stats::rbinom(n_sim, 1, family_y$linkinv(eta_t))
    y_c <- stats::rbinom(n_sim, 1, family_y$linkinv(eta_c))
  }

  # P_med = P(Y(x, M(x)) > Y(x, M(x*))) + 0.5 * P(tie)   (Definition 1)
  mean(y_t > y_c) + 0.5 * mean(y_t == y_c)
}

#' Uniforms for common random numbers in `.pmed_core_simple()`
#'
#' Columns: standardization row, M(x), M(x*), Y | M(x), Y | M(x*).
#'
#' @keywords internal
.crn_uniforms <- function(n_sim) {
  matrix(stats::runif(5L * n_sim), ncol = 5L)
}
