#' Invert an LR test: find where excess(t) = LR(t) - crit crosses zero.
#'
#' Steps outward from `est` in increments of `step` until the statistic exceeds
#' the cutoff, then locates the crossing with uniroot (resolution-independent).
#' `domain = c(lo, hi)` bounds the search; an endpoint pinned at a bound means
#' the interval is open there. The step-out loop is sized to span the whole
#' domain (`max_k` steps), so a crossing anywhere in `[lo, hi]` is reached --
#' a fixed iteration cap could otherwise stop short of a wide domain and return
#' a non-root interior point.
#'
#' @return c(lower = ., upper = .)
#' @keywords internal
.mbco_invert <- function(est, excess, domain, step) {
  lo <- domain[1]
  hi <- domain[2]
  # Enough steps to walk from `est` to either bound (plus a margin), so the
  # `outer` clamp always reaches a bound before the loop is exhausted.
  max_k <- as.integer(ceiling((hi - lo) / step)) + 1L
  endpoint <- function(dir) {
    inner <- est
    for (k in seq_len(max_k)) {
      outer <- min(max(est + dir * step * k, lo), hi)
      if (excess(outer) > 0) {
        return(tryCatch(
          stats::uniroot(excess, sort(c(inner, outer)), tol = step / 100)$root,
          error = function(e) NA_real_
        ))
      }
      inner <- outer
      if (outer <= lo || outer >= hi) {
        return(outer)
      }
    }
    outer
  }
  c(lower = endpoint(-1), upper = endpoint(+1))
}

#' Gaussian exact-null union likelihood
#'
#' Maximize over both reduced models, a = 0 and b = 0, with ML variances.
#'
#' @return Maximized joint log-likelihood, or NA_real_ if the fit fails.
#' @keywords internal
.mbco_ll_null <- function(prep) {
  n <- prep$n
  sse_m_full <- .ols_sse(cbind(prep$Dm, X = prep$X), prep$M)
  sse_y_full <- .ols_sse(cbind(prep$Dy, M = prep$M), prep$Y)
  sse_m_a0 <- .ols_sse(prep$Dm, prep$M)
  sse_y_b0 <- .ols_sse(prep$Dy, prep$Y)
  a0 <- .mbco_gll(sse_m_a0, sse_m_a0 / n, n) +
    .mbco_gll(sse_y_full, sse_y_full / n, n)
  b0 <- .mbco_gll(sse_m_full, sse_m_full / n, n) +
    .mbco_gll(sse_y_b0, sse_y_b0 / n, n)
  if (!all(is.finite(c(a0, b0)))) return(NA_real_)
  max(a0, b0)
}

#' Branch likelihood-ratio statistics for H0: a = 0 and H0: b = 0
#'
#' Nested-model LR statistics with ML variances and the linear nuisance terms
#' profiled, matching `.mbco_ll_null()`. When both are below the interval's
#' chi-square(1) cutoff, neither path is individually significant: the data sit
#' near the double null a = b = 0, where the profile LR at p* = 0.5 is the
#' minimum of these two statistics and the chi-square(1) interval is
#' conservative (coverage about 1 - alpha^2, 0.9975 at the 95% level).
#'
#' @return Named numeric vector `c(a0 = , b0 = )`.
#' @keywords internal
.mbco_branch_lr <- function(prep) {
  n <- prep$n
  gll <- function(sse) .mbco_gll(sse, sse / n, n)
  sse_m_full <- .ols_sse(cbind(prep$Dm, X = prep$X), prep$M)
  sse_y_full <- .ols_sse(cbind(prep$Dy, M = prep$M), prep$Y)
  c(a0 = 2 * (gll(sse_m_full) - gll(.ols_sse(prep$Dm, prep$M))),
    b0 = 2 * (gll(sse_y_full) - gll(.ols_sse(prep$Dy, prep$Y))))
}

#' Profile LR statistic for H0: P_med = p (single-mediator Gaussian MBCO)
#'
#' `-2 * (ll_constrained - ll_free)`; the quantity the interval inverts.
#' `NA` when the constrained profile fails. Used by the composite-null
#' calibration harness (`inst/sim/calibration/calib_mbco_null.R`).
#' @keywords internal
.mbco_lr_at <- function(prep, p) {
  a_sign <- if (sign(prep$a_hat) == 0) 1 else sign(prep$a_hat)
  ll0 <- .mbco_ll_constrained_pmed(prep, stats::qnorm(p), a_sign)
  if (is.na(ll0)) NA_real_ else -2 * (ll0 - prep$ll_free)
}

#' Maximized log-likelihood under the constraint `a*b` = ie0
#'
#' The exact null is the union of a = 0 and b = 0. For nonzero targets,
#' pins b = ie0 / a and profiles linear nuisance terms by OLS.
#' @keywords internal
.mbco_ll_constrained_ie <- function(prep, ie0) {
  if (ie0 == 0) return(.mbco_ll_null(prep))
  results <- lapply(c(-1, 1), function(sa) {
    nll <- function(p) {
      a <- sa * exp(p[1])
      b <- ie0 / a
      if (!all(is.finite(c(a, b))) || a == 0) return(1e100)
      sm <- .mbco_sse(prep$quad_m, a)
      sy <- .mbco_sse(prep$quad_y, b)
      ll <- .mbco_gll(sm, sm / prep$n, prep$n) +
        .mbco_gll(sy, sy / prep$n, prep$n)
      if (is.finite(ll)) -ll else 1e100
    }
    starts <- unique(log(c(max(abs(prep$a_hat), 0.1),
                          max(abs(ie0 / max(abs(prep$b_hat), 0.1)), 0.01), 1)))
    fits <- lapply(starts, function(st) .mbco_optim(st, nll))
    valid <- vapply(fits, function(o) !is.null(o) && o$convergence == 0 &&
      is.finite(o$value) && o$value < 1e99 &&
      -nll(o$par) <= prep$ll_free + 1e-6, logical(1))
    if (!any(valid)) return(NA_real_)
    max(vapply(fits[valid], function(o) -nll(o$par), numeric(1)))
  })
  values <- unlist(results)
  if (any(!is.finite(values))) NA_real_ else max(values)
}

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

  # Use prep's own refit (a_hat, b_hat, Vm_hat, Vy_hat) for the point estimate,
  # NOT extract@a_path / extract@b_path. The CI is inverted around prep$ll_free
  # (from the same lm refits), so centering the estimate on prep guarantees the
  # interval brackets it by construction -- immune to any divergence in how
  # extract@a_path / @b_path were populated (e.g. an SEM-sourced extract).
  a <- prep$a_hat
  b <- prep$b_hat
  delta <- prep$delta
  dhat <- a * delta * b / sqrt(2 * (b^2 * prep$Vm_hat + prep$Vy_hat))
  est <- as.numeric(stats::pnorm(dhat))
  ie_est <- a * b

  a_sign <- sign(a)
  if (a_sign == 0) a_sign <- 1
  crit <- stats::qchisq(ci_level, df = 1)
  excess_pmed <- function(ps) {
    ps <- min(max(ps, 1e-5), 1 - 1e-5)
    ll0 <- .mbco_ll_constrained_pmed(prep, stats::qnorm(ps), a_sign)
    if (is.na(ll0)) {
      return(NA_real_)
    } # failed profiling is not statistical rejection
    -2 * (ll0 - prep$ll_free) - crit
  }
  pmed_scan <- .mbco_invert_verified(est, excess_pmed,
    domain = c(1e-4, 1 - 1e-4), step = 0.01, anchors = 0.5)
  ci <- pmed_scan$bounds
  # Composite-null diagnostic (Phase 8): the interval is conservative when
  # neither path is individually significant.
  pmed_scan$branch_lr <- .mbco_branch_lr(prep)
  pmed_scan$near_double_null <- all(is.finite(pmed_scan$branch_lr)) &&
    all(pmed_scan$branch_lr < crit)

  if (is.finite(prep$sd_ie) && prep$sd_ie > 0) {
    excess_ie <- function(ie0) {
      ll0 <- .mbco_ll_constrained_ie(prep, ie0)
      if (is.na(ll0)) {
        return(NA_real_)
      }
      -2 * (ll0 - prep$ll_free) - crit
    }
    ie_step <- prep$sd_ie / 20
    ie_scan <- .mbco_invert_verified(
      ie_est, excess_ie,
      domain = c(ie_est - 12 * prep$sd_ie, ie_est + 12 * prep$sd_ie),
      step = ie_step, anchors = 0
    )
    ie_ci <- ie_scan$bounds
  } else {
    ie_ci <- c(lower = NA_real_, upper = NA_real_)
    ie_scan <- list(status = "failed", verified = FALSE)
  }

  PmedResult(
    estimate = est,
    ci_lower = unname(ci["lower"]),
    ci_upper = unname(ci["upper"]),
    ci_level = ci_level,
    method = "mbco",
    n_boot = NA_integer_,
    boot_estimates = numeric(0),
    ie_estimate = ie_est,
    ie_ci_lower = unname(ie_ci["lower"]),
    ie_ci_upper = unname(ie_ci["upper"]),
    ie_boot_estimates = numeric(0),
    x_ref = x_ref,
    x_value = x_value,
    source_extract = extract,
    # converged tracks the P_med interval only; the IE interval is reported
    # separately and may be NA on a degenerate design.
    converged = pmed_scan$verified,
    diagnostics = list(pmed = pmed_scan, ie = ie_scan)
  )
}

#' Build design matrices and free-fit summaries for the MBCO optimization.
#'
#' `Dm` is the M-equation design (intercept plus mediator covariates) with the
#'   treatment term `a*X` removed; `Dy` is the Y-equation design (intercept,
#'   treatment, and outcome covariates) with the mediator term `b*M` removed.
#'
#' @keywords internal
.mbco_prep <- function(extract, x_ref, x_value) {
  if (!.is_gaussian(extract@family_m) || !.is_gaussian(extract@family_y)) {
    stop("method = \"mbco\" requires a Gaussian outcome and mediator. ",
      "Use method = \"parametric_bootstrap\" or \"nonparametric_bootstrap\" ",
      "for non-Gaussian models.",
      call. = FALSE
    )
  }
  if (x_value == x_ref) {
    stop("method = \"mbco\" requires x_ref != x_value (a non-degenerate ",
      "treatment contrast).",
      call. = FALSE
    )
  }
  data <- extract@data
  if (is.null(data)) {
    stop("method = \"mbco\" needs the source data to refit constrained models, ",
      "but extract@data is NULL.",
      call. = FALSE
    )
  }

  tx <- extract@treatment
  med <- extract@mediator
  out <- extract@outcome
  cov_m <- setdiff(extract@mediator_predictors, c("(Intercept)", tx))
  cov_y <- setdiff(extract@outcome_predictors, c("(Intercept)", tx, med))

  Dm <- stats::model.matrix(
    stats::reformulate(if (length(cov_m)) cov_m else "1"), data
  )
  Dy <- stats::model.matrix(
    stats::reformulate(c(tx, cov_y)), data
  )

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
    quad_m = .mbco_quadratic(Dm, data[[med]], data[[tx]]),
    quad_y = .mbco_quadratic(Dy, data[[out]], data[[med]]),
    delta = x_value - x_ref,
    a_hat = a_hat, b_hat = b_hat,
    Vm_hat = sum(stats::residuals(fm)^2) / stats::nobs(fm),
    Vy_hat = sum(stats::residuals(fy)^2) / stats::nobs(fy),
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

#' Maximized log-likelihood under the constraint P_med = pnorm(qstar).
#'
#' Profiles out intercepts, the direct effect, and covariate coefficients by OLS,
#' with both compatible coefficient signs explored. Positive Vm and Vy/b^2
#' coordinates enforce the original signed constraint; b is profiled by its
#' quadratic stationary equation. Two numerical coordinates remain.
#'
#' @return Maximized joint log-likelihood, or NA_real_ if the fit fails.
#' @keywords internal
.mbco_ll_constrained_pmed <- function(prep, qstar, a_sign) {
  # a_sign is retained for internal compatibility; both signs are now explored.
  .mbco_fit_pmed(prep, qstar)$ll
}
