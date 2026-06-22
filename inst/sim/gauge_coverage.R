#!/usr/bin/env Rscript
# Coverage study for ward_residual()'s gauge residual W and P_med intervals.
#
# Purpose ---------------------------------------------------------------------
# Validate the CI coverage claims made in R/gauge-pmed.R and the gauge-residual
# vignette: the analytic (symmetric Wald) interval is mildly anti-conservative
# for the right-skewed ratios W = R/OE and P_med = IIE/OE near the null, while
# the bootstrap *percentile* interval restores nominal coverage. This script is
# the validation that turns those theory-asserted figures into measured ones.
#
# Closed-form truth -----------------------------------------------------------
# DGP (E[C]=0, all Gaussian noise):
#   C ~ N(0,1);  A ~ Bernoulli(plogis(g0 + gC*C))
#   M = aM*A + aC*C + N(0,1)
#   Y = bA*A + bM*M + tau*A*M + bC*C + N(0,1)
# The interventional corners are theta(a,a') = bA*a + (bM + tau*a)*aM*a', hence
#   OE  = bA + (bM + tau)*aM
#   IDE = bA
#   IIE = bM*aM
#   R   = tau*aM
#   W   = R/OE,   P_med = IIE/OE     (exact, no Monte-Carlo oracle needed)
#
# Usage -----------------------------------------------------------------------
#   Local smoke:   Rscript inst/sim/gauge_coverage.R --smoke
#   Full local:    Rscript inst/sim/gauge_coverage.R --nsim 1000
#   SLURM array:   one grid cell per task via $SLURM_ARRAY_TASK_ID (1-based);
#                  results written to <outdir>/cell_<id>.rds. Submit an array
#                  sized to nrow(grid()) and concatenate the rds files after.
#
# Output ----------------------------------------------------------------------
# A data frame (one row per grid cell x se_method) with empirical coverage of
# the W and P_med intervals and the Monte-Carlo SE of each coverage estimate
# (sqrt(p(1-p)/nsim)); MCSE is REQUIRED for reading any coverage number.

suppressMessages(library(probmed))

## ---- DGP + closed-form truth ----------------------------------------------
.pars <- list(g0 = -0.2, gC = 0.8, aM = 0.6, aC = 0.4,
              bA = 0.5, bM = 0.7, bC = 0.3)

true_WP <- function(tau, p = .pars) {
  OE  <- p$bA + (p$bM + tau) * p$aM
  IIE <- p$bM * p$aM
  R   <- tau * p$aM
  c(W = R / OE, P = IIE / OE, OE = OE)
}

gen <- function(n, tau, seed, p = .pars) {
  set.seed(seed)
  C <- rnorm(n)
  A <- rbinom(n, 1, plogis(p$g0 + p$gC * C))
  M <- p$aM * A + p$aC * C + rnorm(n)
  Y <- p$bA * A + p$bM * M + tau * A * M + p$bC * C + rnorm(n)
  data.frame(A, M, Y, C)
}

## ---- one replicate: did the CIs cover the closed-form truth? ---------------
one_rep <- function(n, tau, se_method, B, seed) {
  tr <- true_WP(tau)
  fit <- tryCatch(
    ward_residual(gen(n, tau, seed), se_method = se_method, B = B,
                  fieller = FALSE, seed = seed),
    error = function(e) NULL)
  if (is.null(fit)) return(c(W = NA, P = NA))
  c(W = as.numeric(fit@W_ci[1] <= tr["W"] && tr["W"] <= fit@W_ci[2]),
    P = as.numeric(fit@p_med_ci[1] <= tr["P"] && tr["P"] <= fit@p_med_ci[2]))
}

## ---- grid -----------------------------------------------------------------
grid <- function() {
  expand.grid(n = c(500L, 1000L, 2000L, 4000L),
              tau = c(0.0, 0.2, 0.8),       # null, weak, strong interaction
              se_method = c("analytic", "bootstrap"),
              stringsAsFactors = FALSE)
}

run_cell <- function(cell, nsim, B, base_seed = 1000L) {
  reps <- vapply(seq_len(nsim),
                 function(s) one_rep(cell$n, cell$tau, cell$se_method, B,
                                     seed = base_seed + s),
                 numeric(2))
  covW <- mean(reps["W", ], na.rm = TRUE)
  covP <- mean(reps["P", ], na.rm = TRUE)
  mcse <- function(p, m) sqrt(p * (1 - p) / m)
  tr <- true_WP(cell$tau)
  data.frame(n = cell$n, tau = cell$tau, se_method = cell$se_method,
             true_W = unname(tr["W"]), true_P = unname(tr["P"]),
             cov_W = covW, mcse_W = mcse(covW, nsim),
             cov_P = covP, mcse_P = mcse(covP, nsim),
             nsim = nsim, B = B)
}

## ---- entry point ----------------------------------------------------------
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  smoke <- "--smoke" %in% args
  getarg <- function(flag, default) {
    i <- match(flag, args); if (is.na(i)) default else as.integer(args[i + 1L])
  }
  nsim <- if (smoke) 30L else getarg("--nsim", 1000L)
  B    <- if (smoke) 100L else getarg("--B", 400L)
  g <- grid()

  task <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
  if (nzchar(task)) {                       # one cell per array task
    id <- as.integer(task)
    outdir <- Sys.getenv("GAUGE_SIM_OUT", "gauge_sim_out")
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    res <- run_cell(g[id, ], nsim, B)
    saveRDS(res, file.path(outdir, sprintf("cell_%03d.rds", id)))
    message("wrote cell ", id)
    return(invisible(res))
  }

  if (smoke) g <- g[g$n == 500L & g$tau %in% c(0.0, 0.8), ]   # 4 cells only
  out <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) run_cell(g[i, ], nsim, B)))
  print(out, row.names = FALSE)
  invisible(out)
}

if (sys.nframe() == 0L) main()
