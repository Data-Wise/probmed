#' Core P_med Computation Dispatcher
#'
#' @keywords internal
.pmed_compute <- function(extract, x_ref, x_value, method, n_boot, ci_level, seed, ...) {
  switch(method,
    parametric_bootstrap = .pmed_parametric_boot(
      extract, x_ref, x_value, n_boot, ci_level, seed, ...
    ),
    nonparametric_bootstrap = .pmed_nonparametric_boot(
      extract, x_ref, x_value, n_boot, ci_level, seed, ...
    ),
    plugin = .pmed_plugin(
      extract, x_ref, x_value, ...
    )
  )
}

#' Plugin Estimator (Point Estimate Only)
#'
#' @keywords internal
.pmed_plugin <- function(extract, x_ref, x_value, n_sim = 10000, ...) {
  # Extract parameters
  a <- extract@a_path
  b <- extract@b_path
  c_prime <- extract@c_prime
  sigma_y <- extract@sigma_y
  sigma_m <- extract@sigma_m

  # Compute P_med
  pmed_est <- .pmed_core_simple(
    a = a,
    b = b,
    c_prime = c_prime,
    x_ref = x_ref,
    x_value = x_value,
    sigma_y = sigma_y,
    sigma_m = sigma_m,
    n_sim = n_sim
  )

  PmedResult(
    estimate = pmed_est,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    ci_level = NA_real_,
    method = "plugin",
    n_boot = NA_integer_,
    x_ref = x_ref,
    x_value = x_value,
    source_extract = extract,
    converged = TRUE
  )
}

#' Simple P_med Core Computation
#'
#' @keywords internal
.pmed_core_simple <- function(a, b, c_prime, x_ref, x_value, sigma_y = 1, sigma_m = 1, n_sim = 10000) {
  # Handle NA sigma values (for non-Gaussian models)
  # For binary outcomes/mediators, residual variance is not applicable
  if (is.na(sigma_m)) sigma_m <- 0
  if (is.na(sigma_y)) sigma_y <- 0

  # Counterfactual means for M
  mu_m_x <- a * x_value # E[M | X = x]
  mu_m_xref <- a * x_ref # E[M | X = x*]

  # Simulate M under treatment (with residual variance if Gaussian)
  if (sigma_m > 0) {
    m_x <- stats::rnorm(n_sim, mean = mu_m_x, sd = sigma_m)
  } else {
    # For non-Gaussian M (e.g., binary), use expected value
    m_x <- rep(mu_m_x, n_sim)
  }

  # Y counterfactuals (with residual variance if Gaussian)
  if (sigma_y > 0) {
    # Y(x*, M_x) = b*M_x + c'*x* + epsilon_y
    y_xref_mx <- b * m_x + c_prime * x_ref + stats::rnorm(n_sim, 0, sigma_y)

    # Y(x, M_x) = b*M_x + c'*x + epsilon_y
    y_x_mx <- b * m_x + c_prime * x_value + stats::rnorm(n_sim, 0, sigma_y)
  } else {
    # For non-Gaussian Y (e.g., binary), use linear predictor
    # Note: For binary Y, this is on the logit scale
    y_xref_mx <- b * m_x + c_prime * x_ref
    y_x_mx <- b * m_x + c_prime * x_value
  }

  # P_med = P(Y(x*, M_x) > Y(x, M_x))
  pmed <- mean(y_xref_mx > y_x_mx)

  return(pmed)
}
