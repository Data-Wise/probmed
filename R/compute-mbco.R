#' Invert an LR test: find where excess(t) = LR(t) - crit crosses zero.
#'
#' Steps outward from `est` in increments of `step` until the statistic exceeds
#' the cutoff, then locates the crossing with uniroot (resolution-independent).
#' `domain = c(lo, hi)` bounds the search; an endpoint pinned at a bound means
#' the interval is open there.
#'
#' @return c(lower = ., upper = .)
#' @keywords internal
.mbco_invert <- function(est, excess, domain, step) {
  lo <- domain[1]; hi <- domain[2]
  endpoint <- function(dir) {
    inner <- est
    for (k in seq_len(200)) {
      outer <- min(max(est + dir * step * k, lo), hi)
      if (excess(outer) > 0) {
        return(tryCatch(
          stats::uniroot(excess, sort(c(inner, outer)), tol = step / 100)$root,
          error = function(e) NA_real_))
      }
      inner <- outer
      if (outer <= lo || outer >= hi) return(outer)
    }
    outer
  }
  c(lower = endpoint(-1), upper = endpoint(+1))
}

#' Maximized log-likelihood under the constraint a*b = ie0.
#'
#' Pins b = ie0 / a; Vy is free (its MLE is SSE_Y / n). Optimizes over
#' (a, log Vm); all linear nuisance terms are profiled by OLS.
#'
#' @return Maximized joint log-likelihood, or NA_real_ if the fit fails.
#' @keywords internal
.mbco_ll_constrained_ie <- function(prep, ie0) {
  n <- prep$n
  if (abs(ie0) < 1e-12) {  # a*b = 0: take b = 0 (Y ~ Dy, no M)
    sse_m <- .ols_sse(cbind(prep$Dm, X = prep$X), prep$M)
    sse_y <- .ols_sse(prep$Dy, prep$Y)
    return(.mbco_gll(sse_m, sse_m / n, n) + .mbco_gll(sse_y, sse_y / n, n))
  }
  nll <- function(p) {
    a <- p[1]
    if (abs(a) < 1e-8) return(1e10)
    Vm <- exp(p[2])
    b  <- ie0 / a
    sse_m <- .ols_sse(prep$Dm, prep$M - a * prep$X)
    sse_y <- .ols_sse(prep$Dy, prep$Y - b * prep$M)
    -(.mbco_gll(sse_m, Vm, n) + .mbco_gll(sse_y, sse_y / n, n))
  }
  # b = ie0/a makes the a-profile bimodal away from the point estimate, so we
  # warm-start from BOTH basins: a near the free MLE, and a near ie0/b_hat (the
  # basin where b stays near its free value). Keep the better optimum.
  starts <- list(c(prep$a_hat, log(prep$Vm_hat)))
  if (abs(prep$b_hat) > 1e-8) {
    starts[[length(starts) + 1L]] <- c(ie0 / prep$b_hat, log(prep$Vm_hat))
  }
  best <- NA_real_
  for (s in starts) {
    if (!all(is.finite(s))) next
    o <- try(stats::optim(s, nll, method = "Nelder-Mead",
                          control = list(maxit = 3000, reltol = 1e-10)),
             silent = TRUE)
    if (!inherits(o, "try-error") && o$value < 1e9) {
      best <- max(best, -o$value, na.rm = TRUE)
    }
  }
  if (!is.finite(best)) return(NA_real_)
  best
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

  a <- extract@a_path
  b <- extract@b_path
  delta <- prep$delta
  dhat <- a * delta * b / sqrt(2 * (b^2 * prep$Vm_hat + prep$Vy_hat))
  est  <- as.numeric(stats::pnorm(dhat))
  ie_est <- a * b

  a_sign <- sign(a); if (a_sign == 0) a_sign <- 1
  crit <- stats::qchisq(ci_level, df = 1)
  excess_pmed <- function(ps) {
    ps  <- min(max(ps, 1e-5), 1 - 1e-5)
    ll0 <- .mbco_ll_constrained_pmed(prep, stats::qnorm(ps), a_sign)
    if (is.na(ll0)) return(crit + 1e3)          # infeasible -> rejected
    -2 * (ll0 - prep$ll_free) - crit
  }
  ci <- .mbco_invert(est, excess_pmed, domain = c(1e-4, 1 - 1e-4), step = 0.01)

  if (is.finite(prep$sd_ie) && prep$sd_ie > 0) {
    excess_ie <- function(ie0) {
      ll0 <- .mbco_ll_constrained_ie(prep, ie0)
      if (is.na(ll0)) return(crit + 1e3)
      -2 * (ll0 - prep$ll_free) - crit
    }
    ie_step <- prep$sd_ie / 20
    ie_ci <- .mbco_invert(
      ie_est, excess_ie,
      domain = c(ie_est - 12 * prep$sd_ie, ie_est + 12 * prep$sd_ie),
      step = ie_step
    )
  } else {
    ie_ci <- c(lower = NA_real_, upper = NA_real_)
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
    converged = !is.na(ci["lower"]) && !is.na(ci["upper"])
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
  if (x_value == x_ref) {
    stop("method = \"mbco\" requires x_ref != x_value (a non-degenerate ",
         "treatment contrast).", call. = FALSE)
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

#' Maximized log-likelihood under the constraint P_med = pnorm(qstar).
#'
#' Profiles out intercepts, the direct effect, and covariate coefficients by OLS,
#' so only (a, log Vm, log|b|) are optimized numerically (b is sign-pinned to the
#' target branch; Vy is fixed by the P_med constraint). Nelder-Mead is used because
#' the Vy <= 0 infeasible region is a hard penalty wall that defeats gradient methods.
#'
#' @return Maximized joint log-likelihood, or NA_real_ if the fit fails.
#' @keywords internal
.mbco_ll_constrained_pmed <- function(prep, qstar, a_sign) {
  n <- prep$n

  # p* = 0.5: indirect effect is zero (b = 0). Y reduces to its design Dy.
  if (abs(qstar) < 1e-8) {
    sse_m <- .ols_sse(cbind(prep$Dm, X = prep$X), prep$M)  # full M eq: [1, C, X]
    sse_y <- .ols_sse(prep$Dy, prep$Y)                     # Y ~ [1, X, C], no M
    return(.mbco_gll(sse_m, sse_m / n, n) + .mbco_gll(sse_y, sse_y / n, n))
  }

  delta <- prep$delta
  # Constraint requires sign(a*delta*b) = sign(qstar), so sign(a*b) = sign(qstar)*sign(delta).
  sb <- sign(qstar) * a_sign * sign(delta)

  # Optimize over (a, log Vm, log|b|). |b| is a search coordinate, NOT a profiled
  # quantity: b enters both the constraint (which fixes Vy) and the Y residual, so
  # unlike the intercepts / cp / covariate coefficients it cannot be OLS-profiled.
  start <- c(prep$a_hat, log(prep$Vm_hat), log(abs(prep$b_hat) + 1e-3))
  nll <- function(p) {
    a  <- p[1]
    Vm <- exp(p[2])
    b  <- sb * exp(p[3])
    Vy <- (a * delta * b)^2 / (2 * qstar^2) - b^2 * Vm
    # 1e10 penalty wall (infeasible Vy) is kept well above any real NLL so the
    # o$value >= 1e9 check below reliably detects a failed/penalized fit.
    if (!is.finite(Vy) || Vy <= 1e-8) return(1e10)
    sse_m <- .ols_sse(prep$Dm, prep$M - a * prep$X)
    sse_y <- .ols_sse(prep$Dy, prep$Y - b * prep$M)
    -(.mbco_gll(sse_m, Vm, n) + .mbco_gll(sse_y, Vy, n))
  }

  o <- try(stats::optim(start, nll, method = "Nelder-Mead",
                        control = list(maxit = 3000, reltol = 1e-10)),
           silent = TRUE)
  if (inherits(o, "try-error") || o$value >= 1e9) return(NA_real_)
  -o$value
}
