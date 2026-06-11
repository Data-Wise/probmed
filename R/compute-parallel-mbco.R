#' Build design matrices and free-fit summaries for the parallel MBCO fit.
#'
#' Generalizes `.mbco_prep` to k parallel mediators. `Dm_list[[j]]` is the
#' design for mediator j (intercept + its covariates) with the treatment
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
    # Precomputed QR of the FIXED design matrices: the constrained optimizers
    # only change the response (Y - M b, M_j - a_j X), never the design, so
    # decomposing once and reusing qr.resid avoids re-factorizing every nll call.
    qr_Dy = qr(Dy), qr_Dm_list = lapply(Dm_list, qr),
    delta = x_value - x_ref,
    a_hat = a_hat, b_hat = b_hat, Vm_hat = Vm_hat, Vy_hat = Vy_hat,
    ll_free = ll_free, sd_ie = sd_ie
  )
}

#' Residual sum of squares against a precomputed QR decomposition.
#' @keywords internal
.qr_sse <- function(qr_obj, y) {
  sum(qr.resid(qr_obj, y)^2) # qr.resid is base, not stats
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
    sse_mj <- .qr_sse(prep$qr_Dm_list[[j]], prep$M[, j] - a[j] * prep$X)
    ll <- ll + .mbco_gll(sse_mj, Vm[j], n)
  }
  sse_y <- .qr_sse(prep$qr_Dy, prep$Y - as.numeric(prep$M %*% b))
  ll + .mbco_gll(sse_y, Vy, n)
}

#' Deterministic perturbed warm starts for the parallel MBCO optimizers.
#'
#' Returns `base` followed by `n_extra` jittered copies. The jitter draws from a
#' fixed local RNG seed whose prior global state is restored on exit, so the
#' starts (and hence the MBCO interval) stay reproducible and seed-free for the
#' caller. The constrained surfaces are multimodal / penalty-walled, so a single
#' warm start can stick below the true constrained MLE; multistart guards that.
#'
#' @keywords internal
.mbco_parallel_starts <- function(base, n_extra, sd = 0.25) {
  starts <- list(base)
  if (n_extra <= 0L) {
    return(starts)
  }
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(".Random.seed", envir = .GlobalEnv)
      },
      add = TRUE
    )
  }
  set.seed(20260611L)
  for (i in seq_len(n_extra)) {
    starts[[length(starts) + 1L]] <- base + stats::rnorm(length(base), 0, sd)
  }
  starts
}

#' Run Nelder-Mead from several starts, returning the best maximized
#' log-likelihood (`-min value`), or NA if every start fails / is penalized.
#'
#' The first start (the analytic warm start) gets the full iteration budget; the
#' remaining (perturbed) starts are capped low -- they exist only to escape a
#' nearby false basin, so a start that does not converge quickly is abandoned
#' cheaply rather than burning the full budget against the penalty wall.
#'
#' @keywords internal
.mbco_parallel_optim_best <- function(nll, starts, maxit = 1500L,
                                      maxit_extra = 400L) {
  best <- NA_real_
  for (i in seq_along(starts)) {
    s <- starts[[i]]
    if (!all(is.finite(s))) next
    o <- try(
      stats::optim(s, nll,
        method = "Nelder-Mead",
        control = list(
          # reltol 1e-8 is ample for a CI (P_med to ~1e-4) and avoids many extra
          # iterations chasing negligible log-likelihood gains near the optimum.
          maxit = if (i == 1L) maxit else maxit_extra, reltol = 1e-8
        )
      ),
      silent = TRUE
    )
    if (!inherits(o, "try-error") && o$value < 1e9) {
      best <- max(best, -o$value, na.rm = TRUE)
    }
  }
  if (!is.finite(best)) NA_real_ else best
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

  # p* = 0.5  <=>  S = sum(a*b) = 0. For k = 1 this forces b = 0, but for k >= 2
  # the mediators can cancel (e.g. a1*b1 = -a2*b2) with all b_j != 0, which fits
  # the outcome strictly better -- a naive b = 0 fit understates the constrained
  # log-likelihood and inflates the LR statistic near P_med = 0.5. The constraint
  # S = 0 is exactly the ie0 = 0 case of the IE-constrained fit (b_k = (0 - ...) /
  # a_k collapses to b = 0 at k = 1), so reuse it.
  if (abs(qstar) < 1e-8) {
    return(.mbco_ll_constrained_ie_parallel(prep, 0))
  }

  s_target <- sign(qstar) * sign(delta) # required sign of S = sum(a*b)
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

  # 3k coupled coordinates behind a Vy penalty wall: a single start can stick
  # below the constrained MLE (inflating the LR statistic). Multistart.
  base <- c(prep$a_hat, log(prep$Vm_hat), prep$b_hat)
  .mbco_parallel_optim_best(nll, .mbco_parallel_starts(base, n_extra = 3L))
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
    sse_y <- .qr_sse(prep$qr_Dy, prep$Y - as.numeric(prep$M %*% b))
    -.mbco_parallel_ll(prep, a, Vm, b, sse_y / n)
  }

  base <- c(prep$a_hat, log(prep$Vm_hat), prep$b_hat[-k])
  # The b_k = (ie0 - ...)/a_k substitution makes the profile multimodal in each
  # (a_j, b_j) pair. Start from the free MLE, the analytic second basin (a_k
  # tilted to satisfy the constraint), and k-scaled jittered starts to cover the
  # additional basins that appear as k grows.
  starts <- .mbco_parallel_starts(base, n_extra = 2L)
  if (abs(prep$b_hat[k]) > 1e-8) {
    a_alt <- prep$a_hat
    a_alt[k] <- (ie0 - sum(prep$a_hat[-k] * prep$b_hat[-k])) / prep$b_hat[k]
    starts[[length(starts) + 1L]] <- c(a_alt, log(prep$Vm_hat), prep$b_hat[-k])
  }

  .mbco_parallel_optim_best(nll, starts)
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
    # Coarser than the single-mediator step: it only needs to bracket the
    # crossing (uniroot refines to step/100), and each constrained fit here is
    # a multistart optim, so fewer outward steps is a large saving.
    ie_step <- prep$sd_ie / 6
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
