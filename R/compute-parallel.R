#' Closed-form joint P_med for parallel mediators (Gaussian).
#'
#' For k parallel mediators with treatment fixed at the treated level for both
#' potential outcomes, `Y_t - Y_c ~ N(delta * sum(a*b), 2*sum(b^2*Vm) + 2*Vy)`,
#' so the joint P_med is `Phi(delta*sum(a*b) / sqrt(2*sum(b^2*Vm) + 2*Vy))`.
#' Recovers the single-mediator formula at k = 1.
#'
#' @param a_vec,b_vec Length-k path vectors (X -> Mj, Mj -> Y).
#' @param sigma_m_vec Length-k mediator residual SDs.
#' @param sigma_y Outcome residual SD (scalar).
#' @param delta Treatment contrast `x_value - x_ref`.
#' @return Scalar joint P_med in `[0, 1]`.
#' @keywords internal
.pmed_parallel_closed <- function(a_vec, b_vec, sigma_m_vec, sigma_y, delta) {
  num <- delta * sum(a_vec * b_vec)
  den <- sqrt(2 * sum(b_vec^2 * sigma_m_vec^2) + 2 * sigma_y^2)
  as.numeric(stats::pnorm(num / den))
}

#' Require a Gaussian parallel extract (finite residual SDs present).
#'
#' `ParallelMediationData` carries no `family_*` slots, so the Gaussian
#' indicator is "all `sigma_mediators` and `sigma_y` present and finite". A
#' non-Gaussian fit leaves these NA/empty.
#'
#' @keywords internal
.pmed_parallel_require_gaussian <- function(extract) {
  sm <- extract@sigma_mediators
  sy <- extract@sigma_y
  ok <- !is.null(sm) && length(sm) > 0 && all(is.finite(sm)) &&
    !is.null(sy) && length(sy) > 0 && all(is.finite(sy))
  if (!ok) {
    stop("Parallel P_med currently supports a Gaussian outcome and mediators ",
      "only. This `ParallelMediationData` lacks finite residual SDs ",
      "(`sigma_mediators` / `sigma_y`), which indicates a non-Gaussian fit. ",
      "Non-Gaussian parallel mediation is blocked pending family slots in ",
      "medfit's ParallelMediationData.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Dispatch parallel (joint) P_med by inference method.
#'
#' @keywords internal
.pmed_compute_parallel <- function(extract, x_ref, x_value, method,
                                   n_boot, ci_level, seed, ...) {
  switch(method,
    plugin = .pmed_parallel_plugin(
      extract, x_ref, x_value, ...
    ),
    parametric_bootstrap = .pmed_parallel_parametric_boot(
      extract, x_ref, x_value, n_boot, ci_level, seed, ...
    ),
    nonparametric_bootstrap = .pmed_parallel_nonparametric_boot(
      extract, x_ref, x_value, n_boot, ci_level, seed, ...
    ),
    mbco = .pmed_parallel_mbco(
      extract, x_ref, x_value, ci_level = ci_level, ...
    )
  )
}

#' Plugin (point estimate only) for parallel joint P_med.
#'
#' @keywords internal
.pmed_parallel_plugin <- function(extract, x_ref, x_value, ...) {
  .pmed_parallel_require_gaussian(extract)
  a_vec <- extract@a_paths
  b_vec <- extract@b_paths
  delta <- x_value - x_ref
  est <- .pmed_parallel_closed(
    a_vec, b_vec, extract@sigma_mediators, extract@sigma_y, delta
  )
  ie_est <- sum(a_vec * b_vec)

  PmedResult(
    estimate = est,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    ci_level = NA_real_,
    method = "plugin",
    n_boot = NA_integer_,
    boot_estimates = numeric(0),
    ie_estimate = ie_est,
    ie_ci_lower = NA_real_,
    ie_ci_upper = NA_real_,
    ie_boot_estimates = numeric(0),
    x_ref = x_ref,
    x_value = x_value,
    source_extract = extract,
    converged = TRUE
  )
}

#' Names of the structural coefficients (a_j, b_j) in the parallel estimates
#' vector / vcov: `m{j}_<treatment>` for the a-paths, `y_<mediator_j>` for the
#' b-paths (see medfit's parallel extractor).
#'
#' @keywords internal
.pmed_parallel_coef_names <- function(extract) {
  k <- length(extract@mediators)
  list(
    a = paste0("m", seq_len(k), "_", extract@treatment),
    b = paste0("y_", extract@mediators)
  )
}

#' Parametric bootstrap for parallel joint P_med.
#'
#' Draws the structural coefficient block `(a_1..a_k, b_1..b_k)` from
#' `N(estimates, vcov)` (residual SDs held fixed, as in the single-mediator
#' parametric bootstrap) and evaluates the closed-form kernel per draw.
#'
#' @keywords internal
.pmed_parallel_parametric_boot <- function(extract, x_ref, x_value,
                                           n_boot, ci_level, seed, ...) {
  .pmed_parallel_require_gaussian(extract)
  if (!is.null(seed)) set.seed(seed)
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' needed for parametric bootstrap. ",
      "Install with: install.packages('MASS')",
      call. = FALSE
    )
  }

  nm <- .pmed_parallel_coef_names(extract)
  all_nm <- c(nm$a, nm$b)
  mu <- extract@estimates[all_nm]
  Sigma <- extract@vcov[all_nm, all_nm, drop = FALSE]

  theta <- MASS::mvrnorm(n = n_boot, mu = mu, Sigma = Sigma)
  if (n_boot == 1L) {
    theta <- matrix(theta, nrow = 1L, dimnames = list(NULL, all_nm))
  }

  sm <- extract@sigma_mediators
  sy <- extract@sigma_y
  delta <- x_value - x_ref

  boot_estimates <- numeric(n_boot)
  ie_boot_estimates <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    a_vec <- theta[i, nm$a]
    b_vec <- theta[i, nm$b]
    boot_estimates[i] <- .pmed_parallel_closed(a_vec, b_vec, sm, sy, delta)
    ie_boot_estimates[i] <- sum(a_vec * b_vec)
  }

  .pmed_parallel_boot_result(
    boot_estimates, ie_boot_estimates, ci_level, n_boot,
    "parametric_bootstrap", x_ref, x_value, extract
  )
}

#' Nonparametric bootstrap for parallel joint P_med.
#'
#' Resamples rows, refits the k mediator models and the outcome model, and
#' evaluates the closed-form kernel per resample.
#'
#' @keywords internal
.pmed_parallel_nonparametric_boot <- function(extract, x_ref, x_value,
                                              n_boot, ci_level, seed, ...) {
  .pmed_parallel_require_gaussian(extract)
  if (!is.null(seed)) set.seed(seed)

  data <- extract@data
  if (is.null(data)) {
    stop("Data is required for nonparametric bootstrap.", call. = FALSE)
  }
  n <- nrow(data)
  tx <- extract@treatment
  meds <- extract@mediators
  out <- extract@outcome
  k <- length(meds)

  med_forms <- lapply(seq_len(k), function(j) {
    preds <- setdiff(extract@mediator_predictors[[j]], "(Intercept)")
    if (length(preds) == 0L) preds <- "1"
    stats::reformulate(preds, response = meds[j])
  })
  preds_y <- setdiff(extract@outcome_predictors, "(Intercept)")
  if (length(preds_y) == 0L) preds_y <- "1"
  form_y <- stats::reformulate(preds_y, response = out)

  delta <- x_value - x_ref
  boot_estimates <- numeric(n_boot)
  ie_boot_estimates <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    bd <- data[idx, , drop = FALSE]
    fits_m <- lapply(med_forms, function(f) stats::lm(f, bd))
    fy <- stats::lm(form_y, bd)
    cy <- stats::coef(fy)
    a_vec <- vapply(fits_m, function(fm) unname(stats::coef(fm)[tx]), numeric(1))
    b_vec <- vapply(meds, function(m) unname(cy[m]), numeric(1))
    sm_vec <- vapply(fits_m, function(fm) stats::sigma(fm), numeric(1))
    sy <- stats::sigma(fy)
    boot_estimates[i] <- .pmed_parallel_closed(a_vec, b_vec, sm_vec, sy, delta)
    ie_boot_estimates[i] <- sum(a_vec * b_vec)
  }

  .pmed_parallel_boot_result(
    boot_estimates, ie_boot_estimates, ci_level, n_boot,
    "nonparametric_bootstrap", x_ref, x_value, extract
  )
}

#' Assemble a bootstrap `PmedResult` (shared by both parallel bootstraps).
#'
#' @keywords internal
.pmed_parallel_boot_result <- function(boot_estimates, ie_boot_estimates,
                                       ci_level, n_boot, method,
                                       x_ref, x_value, extract) {
  alpha <- 1 - ci_level
  probs <- c(alpha / 2, 1 - alpha / 2)
  ci <- stats::quantile(boot_estimates, probs = probs, na.rm = TRUE)
  ie_ci <- stats::quantile(ie_boot_estimates, probs = probs, na.rm = TRUE)

  PmedResult(
    estimate = mean(boot_estimates, na.rm = TRUE),
    ci_lower = unname(ci[1]),
    ci_upper = unname(ci[2]),
    ci_level = ci_level,
    method = method,
    n_boot = as.integer(n_boot),
    boot_estimates = boot_estimates,
    ie_estimate = mean(ie_boot_estimates, na.rm = TRUE),
    ie_ci_lower = unname(ie_ci[1]),
    ie_ci_upper = unname(ie_ci[2]),
    ie_boot_estimates = ie_boot_estimates,
    x_ref = x_ref,
    x_value = x_value,
    source_extract = extract,
    converged = TRUE
  )
}
