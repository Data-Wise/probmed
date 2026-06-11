#' MBCO Confidence Interval for P_med (Gaussian single-mediator model)
#'
#' Implements the Model-Based Constrained Optimization interval (Tofighi &
#' Kelley, 2020) by likelihood-ratio inversion. Gaussian outcome AND mediator
#' only; covariates in either equation are supported as profiled nuisance
#' parameters. Deterministic (no resampling).
#'
#' @keywords internal
.pmed_mbco <- function(extract, x_ref, x_value, ci_level = 0.95, ...) {
  prep <- .mbco_prep(extract, x_ref, x_value)

  a <- extract@a_path
  b <- extract@b_path
  Vm <- prep$Vm_hat
  Vy <- prep$Vy_hat
  delta <- prep$delta

  dhat <- a * delta * b / sqrt(2 * (b^2 * Vm + Vy))
  est  <- as.numeric(stats::pnorm(dhat))
  ie_est <- a * b

  PmedResult(
    estimate = est,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    ci_level = NA_real_,
    method = "mbco",
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

#' Build design matrices and free-fit summaries for the MBCO optimization.
#'
#' Dm = [intercept, mediator covariates] (the M-equation design with a*X removed).
#' Dy = [intercept, treatment, outcome covariates] (the Y-equation design with b*M removed).
#'
#' @keywords internal
.mbco_prep <- function(extract, x_ref, x_value) {
  if (!.is_gaussian(extract@family_m) || !.is_gaussian(extract@family_y)) {
    stop("method = \"mbco\" requires a Gaussian outcome and mediator. ",
         "Use method = \"parametric_bootstrap\" or \"nonparametric_bootstrap\" ",
         "for non-Gaussian models.", call. = FALSE)
  }
  data <- extract@data
  if (is.null(data)) {
    stop("method = \"mbco\" needs the source data to refit constrained models, ",
         "but extract@data is NULL.", call. = FALSE)
  }

  tx  <- extract@treatment
  med <- extract@mediator
  out <- extract@outcome
  cov_m <- setdiff(extract@mediator_predictors, c("(Intercept)", tx))
  cov_y <- setdiff(extract@outcome_predictors, c("(Intercept)", tx, med))

  Dm <- stats::model.matrix(
    stats::reformulate(if (length(cov_m)) cov_m else "1"), data)
  Dy <- stats::model.matrix(
    stats::reformulate(c(tx, cov_y)), data)

  # Free fits (full M and Y equations) for warm starts, point estimates, scale.
  fm <- stats::lm(stats::reformulate(c(tx, cov_m), med), data)
  fy <- stats::lm(stats::reformulate(c(tx, med, cov_y), out), data)
  a_hat <- unname(stats::coef(fm)[tx])
  b_hat <- unname(stats::coef(fy)[med])
  sm <- summary(fm)$coefficients
  sy <- summary(fy)$coefficients
  se_a <- sm[tx, "Std. Error"]
  se_b <- sy[med, "Std. Error"]

  list(
    M = data[[med]], Y = data[[out]], X = data[[tx]],
    Dm = Dm, Dy = Dy, n = nrow(data),
    delta = x_value - x_ref,
    a_hat = a_hat, b_hat = b_hat,
    Vm_hat = stats::sigma(fm)^2, Vy_hat = stats::sigma(fy)^2,
    ll_free = as.numeric(stats::logLik(fm)) + as.numeric(stats::logLik(fy)),
    # delta-method SD of a*b, for the IE inversion search scale
    sd_ie = sqrt(b_hat^2 * se_a^2 + a_hat^2 * se_b^2)
  )
}

#' Residual sum of squares of an OLS fit of y on design matrix D.
#' @keywords internal
.ols_sse <- function(D, y) {
  fit <- stats::lm.fit(D, y)
  sum(fit$residuals^2)
}

#' Gaussian profile log-likelihood given SSE, variance V, and n rows.
#' @keywords internal
.mbco_gll <- function(sse, V, n) {
  -n / 2 * log(2 * pi * V) - sse / (2 * V)
}
