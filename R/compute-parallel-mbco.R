#' Build design matrices and free-fit summaries for the parallel MBCO fit.
#'
#' Generalizes `.mbco_prep` to k parallel mediators. `Dm_list[[j]]` is the j-th
#' mediator design (intercept + that mediator's covariates) with the treatment
#' term removed; `Dy` is the outcome design (intercept + treatment + outcome
#' covariates) with all mediator terms removed. The b-paths multiply the columns
#' of `M` (n x k).
#'
#' @keywords internal
.mbco_prep_parallel <- function(extract, x_ref, x_value) {
  .pmed_parallel_require_gaussian(extract)
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
  meds <- extract@mediators
  out <- extract@outcome
  k <- length(meds)
  cov_y <- setdiff(extract@outcome_predictors, c("(Intercept)", tx, meds))

  cov_m_list <- lapply(seq_len(k), function(j) {
    setdiff(extract@mediator_predictors[[j]], c("(Intercept)", tx))
  })
  Dm_list <- lapply(cov_m_list, function(cm) {
    stats::model.matrix(stats::reformulate(if (length(cm)) cm else "1"), data)
  })
  Dy <- stats::model.matrix(stats::reformulate(c(tx, cov_y)), data)

  fm_list <- lapply(seq_len(k), function(j) {
    stats::lm(stats::reformulate(c(tx, cov_m_list[[j]]), meds[j]), data)
  })
  fy <- stats::lm(stats::reformulate(c(tx, meds, cov_y), out), data)

  cy <- stats::coef(fy)
  a_hat <- vapply(fm_list, function(fm) unname(stats::coef(fm)[tx]), numeric(1))
  b_hat <- vapply(meds, function(m) unname(cy[m]), numeric(1))
  Vm_hat <- vapply(fm_list, function(fm) stats::sigma(fm)^2, numeric(1))
  Vy_hat <- stats::sigma(fy)^2
  ll_free <- sum(vapply(
    fm_list, function(fm) as.numeric(stats::logLik(fm)), numeric(1)
  )) + as.numeric(stats::logLik(fy))

  # Delta-method SD of S = sum(a_j b_j), for the IE inversion search scale only
  # (independence approximation; the scale need not be exact).
  se_a <- vapply(seq_len(k), function(j) {
    summary(fm_list[[j]])$coefficients[tx, "Std. Error"]
  }, numeric(1))
  sy_co <- summary(fy)$coefficients
  se_b <- vapply(meds, function(m) sy_co[m, "Std. Error"], numeric(1))
  sd_ie <- sqrt(sum(b_hat^2 * se_a^2 + a_hat^2 * se_b^2))

  list(
    M = as.matrix(data[, meds, drop = FALSE]),
    Y = data[[out]], X = data[[tx]],
    Dm_list = Dm_list, Dy = Dy, n = nrow(data), k = k,
    delta = x_value - x_ref,
    a_hat = a_hat, b_hat = b_hat, Vm_hat = Vm_hat, Vy_hat = Vy_hat,
    ll_free = ll_free, sd_ie = sd_ie
  )
}

#' Joint Gaussian log-likelihood for the parallel free model at given params.
#'
#' Each mediator residual `M[,j] - a_j X` is OLS-profiled on `Dm_list[[j]]`;
#' the outcome residual `Y - M b` is OLS-profiled on `Dy`. Returns the summed
#' maximized log-likelihood for the supplied `(a, Vm, b, Vy)`.
#'
#' @keywords internal
.mbco_parallel_ll <- function(prep, a, Vm, b, Vy) {
  n <- prep$n
  ll <- 0
  for (j in seq_len(prep$k)) {
    sse_mj <- .ols_sse(prep$Dm_list[[j]], prep$M[, j] - a[j] * prep$X)
    ll <- ll + .mbco_gll(sse_mj, Vm[j], n)
  }
  sse_y <- .ols_sse(prep$Dy, prep$Y - as.numeric(prep$M %*% b))
  ll + .mbco_gll(sse_y, Vy, n)
}

#' Maximized log-likelihood under the constraint joint P_med = pnorm(qstar).
#'
#' The constraint fixes `Vy = (delta*S)^2 / (2 qstar^2) - sum(b^2 Vm)` with
#' `S = sum(a*b)`, and requires `sign(delta*S) = sign(qstar)`. Optimizes over
#' `(a_1..a_k, log Vm_1..k, b_1..k)` (3k coords) with intercepts / direct effect
#' / covariates OLS-profiled. Nelder-Mead (infeasible Vy is a penalty wall).
#'
#' @keywords internal
.mbco_ll_constrained_pmed_parallel <- function(prep, qstar) {
  n <- prep$n
  k <- prep$k
  delta <- prep$delta

  # p* = 0.5: indirect effect S = 0; take b = 0 (Y ~ Dy, no mediators), each
  # mediator equation full (intercept + treatment + covariates).
  if (abs(qstar) < 1e-8) {
    sse_y <- .ols_sse(prep$Dy, prep$Y)
    ll <- .mbco_gll(sse_y, sse_y / n, n)
    for (j in seq_len(k)) {
      sse_mj <- .ols_sse(cbind(prep$Dm_list[[j]], X = prep$X), prep$M[, j])
      ll <- ll + .mbco_gll(sse_mj, sse_mj / n, n)
    }
    return(ll)
  }

  s_target <- sign(qstar) * sign(delta) # required sign of S = sum(a*b)
  start <- c(prep$a_hat, log(prep$Vm_hat), prep$b_hat)
  nll <- function(p) {
    a <- p[seq_len(k)]
    Vm <- exp(p[k + seq_len(k)])
    b <- p[2 * k + seq_len(k)]
    S <- sum(a * b)
    if (abs(S) > 1e-10 && sign(S) != s_target) {
      return(1e10)
    }
    Vy <- (delta * S)^2 / (2 * qstar^2) - sum(b^2 * Vm)
    if (!is.finite(Vy) || Vy <= 1e-8) {
      return(1e10)
    }
    -.mbco_parallel_ll(prep, a, Vm, b, Vy)
  }

  o <- try(
    stats::optim(start, nll,
      method = "Nelder-Mead",
      control = list(maxit = 5000, reltol = 1e-10)
    ),
    silent = TRUE
  )
  if (inherits(o, "try-error") || o$value >= 1e9) {
    return(NA_real_)
  }
  -o$value
}

#' Maximized log-likelihood under the constraint sum(a*b) = ie0 (parallel).
#'
#' Pins the last b-path via `b_k = (ie0 - sum_{j<k} a_j b_j) / a_k`; Vy is free
#' (MLE = SSE_Y / n). Optimizes over `(a_1..k, log Vm_1..k, b_1..k-1)`. Because
#' the b_k = .../a_k substitution makes the profile multimodal away from the
#' point estimate, warm-starts from two basins (as the single-mediator IE fit).
#'
#' @keywords internal
.mbco_ll_constrained_ie_parallel <- function(prep, ie0) {
  n <- prep$n
  k <- prep$k

  nll <- function(p) {
    a <- p[seq_len(k)]
    Vm <- exp(p[k + seq_len(k)])
    b_head <- if (k > 1L) p[2 * k + seq_len(k - 1L)] else numeric(0)
    if (abs(a[k]) < 1e-8) {
      return(1e10)
    }
    b_k <- (ie0 - sum(a[-k] * b_head)) / a[k]
    b <- c(b_head, b_k)
    sse_y <- .ols_sse(prep$Dy, prep$Y - as.numeric(prep$M %*% b))
    -.mbco_parallel_ll(prep, a, Vm, b, sse_y / n)
  }

  starts <- list(c(prep$a_hat, log(prep$Vm_hat), prep$b_hat[-k]))
  # Second basin: a_k near (ie0 - sum_{j<k} a_j b_j)/b_k, keeping b near free.
  if (abs(prep$b_hat[k]) > 1e-8) {
    a_alt <- prep$a_hat
    a_alt[k] <- (ie0 - sum(prep$a_hat[-k] * prep$b_hat[-k])) / prep$b_hat[k]
    starts[[length(starts) + 1L]] <- c(a_alt, log(prep$Vm_hat), prep$b_hat[-k])
  }

  best <- NA_real_
  for (s in starts) {
    if (!all(is.finite(s))) next
    o <- try(
      stats::optim(s, nll,
        method = "Nelder-Mead",
        control = list(maxit = 5000, reltol = 1e-10)
      ),
      silent = TRUE
    )
    if (!inherits(o, "try-error") && o$value < 1e9) {
      best <- max(best, -o$value, na.rm = TRUE)
    }
  }
  if (!is.finite(best)) {
    return(NA_real_)
  }
  best
}

#' MBCO confidence interval for parallel joint P_med (Gaussian).
#'
#' Deterministic likelihood-ratio inversion (Tofighi & Kelley, 2020) for the
#' joint parallel estimand and for the total indirect effect `sum(a*b)`, reusing
#' the generic `.mbco_invert` engine. Gaussian outcome and mediators only.
#'
#' @keywords internal
.pmed_parallel_mbco <- function(extract, x_ref, x_value, ci_level = 0.95, ...) {
  prep <- .mbco_prep_parallel(extract, x_ref, x_value)

  a <- prep$a_hat
  b <- prep$b_hat
  delta <- prep$delta
  s_hat <- sum(a * b)
  dhat <- delta * s_hat / sqrt(2 * sum(b^2 * prep$Vm_hat) + 2 * prep$Vy_hat)
  est <- as.numeric(stats::pnorm(dhat))
  ie_est <- s_hat

  crit <- stats::qchisq(ci_level, df = 1)
  excess_pmed <- function(ps) {
    ps <- min(max(ps, 1e-5), 1 - 1e-5)
    ll0 <- .mbco_ll_constrained_pmed_parallel(prep, stats::qnorm(ps))
    if (is.na(ll0)) {
      return(crit + 1e3)
    }
    -2 * (ll0 - prep$ll_free) - crit
  }
  ci <- .mbco_invert(est, excess_pmed, domain = c(1e-4, 1 - 1e-4), step = 0.01)

  if (is.finite(prep$sd_ie) && prep$sd_ie > 0) {
    excess_ie <- function(ie0) {
      ll0 <- .mbco_ll_constrained_ie_parallel(prep, ie0)
      if (is.na(ll0)) {
        return(crit + 1e3)
      }
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
