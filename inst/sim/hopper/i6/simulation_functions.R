# simulation_functions.R -- estimators, bootstrap, and one-replication driver.
# Sourced by run_simulation.R after config.R (and mbco_pmed.R for the MBCO CI).
#
# P_med is estimated by Monte Carlo g-computation -- the same estimator
# probmed::pmed() implements. For speed and zero load-time dependencies the
# simulation uses this self-contained version; the empirical example (memory_exp)
# uses probmed directly. The MBCO interval is the self-contained closed-form
# routine in code/mbco_pmed.R (the same procedure as RMediation::mbco; Tofighi & Kelley, 2020).

# ---- competing point estimators --------------------------------------------
pm_est <- function(data) {                       # proportion mediated a*b/c
  a  <- coef(lm(M ~ X, data))["X"]
  b  <- coef(lm(Y ~ X + M, data))["M"]
  cc <- coef(lm(Y ~ X, data))["X"]
  as.numeric(a * b / cc)
}

beta_std_est <- function(data) {                 # standardized indirect a*b/SD(Y)
  a <- coef(lm(M ~ X, data))["X"]
  b <- coef(lm(Y ~ X + M, data))["M"]
  as.numeric(a * b / stats::sd(data$Y))
}

# ---- P_med point estimate (MC g-computation; treatment fixed at X = 1) ------
# Draws mediator values under control and treatment (independent: cross-world),
# propagates them through the outcome model, and compares by rank.
# `u` (optional, K x 4 uniforms from crn_uniforms()) replaces the fresh draws
# by inverse-CDF transforms: common random numbers across bootstrap resamples.
pmed_point <- function(data, family_y = "gaussian", K = 5000, u = NULL) {
  fm  <- lm(M ~ X, data)
  b0m <- coef(fm)["(Intercept)"]; am <- coef(fm)["X"]; sm <- sigma(fm)
  if (is.null(u)) {
    M0 <- rnorm(K, b0m, sm); M1 <- rnorm(K, b0m + am, sm)
  } else {
    M0 <- b0m + sm * qnorm(u[, 1]); M1 <- b0m + am + sm * qnorm(u[, 2])
  }
  if (family_y == "gaussian") {
    fy <- lm(Y ~ X + M, data); cy <- coef(fy); sy <- sigma(fy)
    lin <- function(m) cy["(Intercept)"] + cy["X"] * 1 + cy["M"] * m
    if (is.null(u)) {
      Y0 <- lin(M0) + rnorm(K, 0, sy); Y1 <- lin(M1) + rnorm(K, 0, sy)
    } else {
      Y0 <- lin(M0) + sy * qnorm(u[, 3]); Y1 <- lin(M1) + sy * qnorm(u[, 4])
    }
    result <- mean(Y1 > Y0)
    attr(result, "sep_warning") <- FALSE
    result
  } else {
    # Finding #21: a binary outcome model with quasi/complete separation (common
    # at small N, or under the rare-binary scenarios) makes glm() report the
    # classic "fitted probabilities numerically 0 or 1 occurred" warning and/or
    # fail to converge; the coefficients are then effectively noise (infinite
    # SE), and pmed_point() previously used them unflagged. Detect both signals
    # per fit rather than silently trusting whatever glm() returns.
    sep <- FALSE
    fy <- withCallingHandlers(
      glm(Y ~ X + M, data, family = binomial),
      warning = function(w) {
        if (grepl("fitted (probabilities|prob\\.?) numerically 0 or 1",
                   conditionMessage(w)))
          sep <<- TRUE
        invokeRestart("muffleWarning")
      })
    if (!isTRUE(fy$converged)) sep <- TRUE
    cy <- coef(fy)
    p  <- function(m) plogis(cy["(Intercept)"] + cy["X"] * 1 + cy["M"] * m)
    if (is.null(u)) {
      Y0 <- rbinom(K, 1, p(M0)); Y1 <- rbinom(K, 1, p(M1))
    } else {
      Y0 <- as.integer(u[, 3] < p(M0)); Y1 <- as.integer(u[, 4] < p(M1))
    }
    result <- mean(Y1 > Y0) + 0.5 * mean(Y1 == Y0)
    attr(result, "sep_warning") <- sep
    result
  }
}

crn_uniforms <- function(K) matrix(runif(4L * K), ncol = 4L)

estimate_all <- function(data, family_y, K, u = NULL) {
  pmed_val <- pmed_point(data, family_y, K, u)
  out <- c(PM       = pm_est(data),
           beta_std = beta_std_est(data),
           pmed     = as.numeric(pmed_val))
  attr(out, "sep_warning") <- isTRUE(attr(pmed_val, "sep_warning"))
  out
}

# ---- one replication: point estimates + percentile-bootstrap CIs -----------
# do_mbco adds an MBCO CI for P_med on continuous-outcome scenarios (slow).
# sep_warnings counts glm() separation/non-convergence events (finding #21)
# across the point estimate plus all n_boot resamples, so a per-cell rate is
# recoverable from the aggregated output without re-running anything.
run_replication <- function(scenario, N, K, n_boot, do_mbco = FALSE) {
  d  <- scenario$gen(N); fy <- scenario$family_y
  est <- estimate_all(d, fy, K)
  sep_count <- as.integer(isTRUE(attr(est, "sep_warning")))
  # Monte Carlo control (MC_CONTROL == "crn", the default; probmed PR #46): every
  # resample evaluates P_med on one shared u, re-centered on the full-sample fit
  # with the same u (a control variate), so the K-pair simulation error neither
  # widens nor shifts the interval. The point estimate uses 100 * K pairs.
  # MC_CONTROL == "fresh" reproduces the original runs (fresh draws per resample).
  crn <- identical(MC_CONTROL, "crn")
  K_POINT <- K
  # K_BOOT pairs per resample evaluation: with shared u, too few pairs make the
  # interval width vary across datasets and undercover at the null (Mech. 6,
  # N = 1000: 0.933 at K = 5,000 vs 0.951 with exact per-resample P_med).
  if (crn) {
    K <- K_BOOT
    u <- crn_uniforms(K)
    anchor <- as.numeric(pmed_point(d, fy, K, u))
    est["pmed"] <- as.numeric(pmed_point(d, fy, 100L * K_POINT))
  } else {
    u <- NULL
  }
  boot <- replicate(n_boot, {
    db <- d[sample.int(N, N, replace = TRUE), ]
    r <- tryCatch(estimate_all(db, fy, K, u), error = function(e) rep(NA_real_, length(est)))
    if (is.numeric(r) && isTRUE(attr(r, "sep_warning"))) sep_count <<- sep_count + 1L
    if (crn) r[3] <- min(max(r[3] - anchor + est[["pmed"]], 0), 1)
    r
  })
  ci  <- apply(boot, 1, quantile, c(.025, .975), na.rm = TRUE)   # 2 x 3
  out <- data.frame(estimator = names(est), estimate = as.numeric(est),
                    ci_lower = ci[1, ], ci_upper = ci[2, ],
                    method = "percentile_boot", sep_warnings = sep_count, row.names = NULL)
  if (do_mbco && fy == "gaussian" && exists("pmed_mbco_ci")) {
    mb <- tryCatch(pmed_mbco_ci(d),
                   error = function(e) NULL)
    if (!is.null(mb))
      out <- rbind(out, data.frame(estimator = "pmed", estimate = mb$estimate,
                   ci_lower = mb$ci_lower, ci_upper = mb$ci_upper,
                   method = "mbco", sep_warnings = sep_count, row.names = NULL))
  }
  out
}

# ---- true P_med per scenario (the DGM estimand, from its potential outcomes) -
# Draws the cross-world potential outcomes Y(1, M(1)) and Y(1, M(0)) directly from
# the true DGM (scenario$po) and forms P_med, so the truth is the estimand itself
# rather than a working-model fit -- exact even where the estimation model is
# misspecified (the s5 mixture; the logistic fit on the probit binary scenarios).
true_pmed <- function(scenario, Nbig = 2e6) {
  po <- scenario$po(Nbig)
  mean(po$Y1 > po$Y0) + 0.5 * mean(po$Y1 == po$Y0)
}
