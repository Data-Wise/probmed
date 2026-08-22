#!/usr/bin/env Rscript
# Decompose the under-dispersion of ward_residual()'s analytic SE for W = R/OE.
#
# Purpose ---------------------------------------------------------------------
# The gauge coverage grid (inst/sim/results/gauge_boot_coverage_nsim2000.csv)
# shows seW_an_ratio = mean(seW_an)/empSD = 0.65-0.83 in ALL 8 manuscript cells:
# the influence-function SE the package reports understates the true sampling
# SD of W-hat by 17-35%. Every interval built from the influence matrix `phi` --
# the Wald arm, any phi-reweighted bootstrap, and the Fieller set -- inherits
# that factor. This script measures WHERE the missing variance comes from.
# Motivation and the refutation that made it the blocking measurement:
# docs/specs/REVIEW-2026-08-22-adversarial-refutation.md.
#
# The identity ----------------------------------------------------------------
# Across datasets, with the as-shipped single random fold partition:
#   V_cf  = Var(W_cf)               true sampling variance of the shipped estimator
#   S_an  = E[seW_an^2]             the analytic variance the package reports
#   gap G = V_cf - S_an = F + N + R, where
#   F = E_data[ Var_partition(W_cf | data) ]    FOLD-SPLIT: same data, re-drawn folds.
#       .corner_fit() draws `folds <- sample(...)` on every call (R/corner.R:29),
#       so the shipped W carries partition noise the IF SE does not see.
#   N = (V_cf - F) - V_or                       NUISANCE ESTIMATION: variance of the
#       partition-averaged cross-fit estimator (law of total variance) minus the
#       variance of the ORACLE estimator that plugs the true nuisances into the
#       same corner EIF. No cross-fitting is needed with known nuisances.
#   R = V_or - S_an                             REMAINDER: the IF/delta-method
#       variance formula itself, at this n, for a ratio. Sub-split as
#       (V_or - S_or) + (S_or - S_an): formula inadequacy with oracle nuisances,
#       plus the effect of plugging estimated nuisances into the formula.
# F + N + R = G by construction; the shares are reported.
#
# Oracle nuisances (closed form from the DGP) ----------------------------------
#   pi(C)      = plogis(-0.2 + 0.8 C)
#   M | A=a,C  ~ N(0.6 s a + 0.4 C, 1)          =>  q(M,C) by Bayes' rule
#   mu(a,M,C)  = lp, or plogis(lp) if binY;  lp = 0.5 s a + 0.7 M + tau s a M + 0.3 C
#   eta(a,a',C)= E[mu(a,M,C) | A=a', C]: linear in M when continuous; a 1-D
#                Gaussian expectation under Gauss-Hermite quadrature when binary.
# With s = 1 and tau in {0, 0.8} this is the manuscript DGP (gauge-pmed.qmd:207)
# exactly; tau = 0.4 with s < 1 is the weak-ID grid's DGP
# (inst/sim/hopper/run_weakid_validation.R:58). Exact truth reuses the weak-ID
# grid's closed form (never Monte Carlo -- see the warning there).
#
# Positive controls (the run is wrong if these fail) ---------------------------
#   * E[phi_oracle[, j]] == theta_exact(a, a')    asserted on a 2e5 sample.
#   * Gauss-Hermite: E[Z^2] == 1, E[plogis(Z)] == 1/2.
#   * Manuscript cells must REPRODUCE seW_an_ratio 0.65-0.83 and Wald coverage
#     0.86-0.91 from the published grid; if they do not, the wiring is wrong.
#
# Usage -----------------------------------------------------------------------
#   Smoke (cells 1,5,9; nrep 20; 4 partitions):  Rscript inst/sim/phi_decomposition.R --smoke
#   One cell locally:                            Rscript inst/sim/phi_decomposition.R --cell 1 --nrep 250
#   SLURM array, one cell per task:              $SLURM_ARRAY_TASK_ID selects the cell
#   Collate all cell files into one table:       Rscript inst/sim/phi_decomposition.R --collate
# Run from the package root. Uses the in-tree package via pkgload when a
# DESCRIPTION is present, else the installed probmed. Per the hopper rule,
# pilot with `sbatch --array=1-1` before a full array.
#
# Output ----------------------------------------------------------------------
#   <outdir>/phi_decomp_cell_<id>.rds : list(summary = 1-row data.frame,
#                                            rows = per-rep data.frame,
#                                            W_cf = nrep x R_part matrix, ...)
#   --collate writes <outdir>/phi_decomp_summary.csv and prints the table.
# outdir = $GAUGE_SIM_OUT/phi_decomp; GAUGE_SIM_OUT defaults to gauge_sim_out,
# the same scratch location the other gauge sim scripts use (PR #26). Promote a
# finished run's summary CSV into inst/sim/results/ deliberately, by hand.

## ---- arguments --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
flag <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
SMOKE   <- isTRUE(flag("smoke", FALSE))
COLLATE <- isTRUE(flag("collate", FALSE))
NREP    <- as.integer(flag("nrep",  if (SMOKE) 20L else 250L))
R_PART  <- as.integer(flag("rpart", if (SMOKE) 4L else 10L))
K       <- 5L
OUTDIR  <- file.path(Sys.getenv("GAUGE_SIM_OUT", "gauge_sim_out"), "phi_decomp")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ---- package: in-tree if possible, else installed ---------------------------
if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
  corner_fit <- .corner_fit
} else {
  suppressMessages(library(probmed))
  corner_fit <- probmed:::.corner_fit
}

## ---- DGP ---------------------------------------------------------------------
gen <- function(n, s, tau, binY) {
  C <- rnorm(n); A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
  M <- 0.6 * s * A + 0.4 * C + rnorm(n)
  lp <- 0.5 * s * A + 0.7 * M + tau * s * A * M + 0.3 * C
  Y <- if (binY) rbinom(n, 1, plogis(lp)) else lp + rnorm(n)
  data.frame(A = A, M = M, Y = Y, C = C)
}

## ---- Gauss-Hermite for E[f(Z)], Z ~ N(0,1) (Golub-Welsch) -------------------
gh_nodes <- function(m) {
  k <- seq_len(m - 1L)
  J <- matrix(0, m, m); J[cbind(k, k + 1L)] <- J[cbind(k + 1L, k)] <- sqrt(k / 2)
  e <- eigen(J, symmetric = TRUE)
  list(x = sqrt(2) * e$values, w = e$vectors[1, ]^2)   # weights sum to 1
}
GH <- gh_nodes(40L)
stopifnot(abs(sum(GH$w * GH$x^2) - 1) < 1e-10,
          abs(sum(GH$w * plogis(GH$x)) - 0.5) < 1e-10)

## ---- oracle nuisances + oracle influence matrix ------------------------------
mu_oracle <- function(a, M, C, s, tau, binY) {
  lp <- 0.5 * s * a + 0.7 * M + tau * s * a * M + 0.3 * C
  if (binY) plogis(lp) else lp
}
eta_oracle <- function(a, ap, C, s, tau, binY) {
  b <- 0.7 + tau * s * a
  base <- 0.5 * s * a + b * (0.6 * s * ap + 0.4 * C) + 0.3 * C   # lp at e_M = 0
  if (!binY) return(base)
  as.vector(plogis(outer(base, b * GH$x, "+")) %*% GH$w)          # E over e_M ~ N(0,1)
}
phi_oracle <- function(d, s, tau, binY) {
  n <- nrow(d); A <- d$A; M <- d$M; Y <- d$Y; C <- d$C
  p1 <- plogis(-0.2 + 0.8 * C)
  f1 <- dnorm(M, 0.6 * s + 0.4 * C, 1); f0 <- dnorm(M, 0.4 * C, 1)
  q1 <- p1 * f1 / (p1 * f1 + (1 - p1) * f0)
  pa <- function(z) if (z == 1) p1 else 1 - p1
  qa <- function(z) if (z == 1) q1 else 1 - q1
  nm <- c("11", "10", "01", "00"); cor <- list(c(1, 1), c(1, 0), c(0, 1), c(0, 0))
  phi <- matrix(0, n, 4, dimnames = list(NULL, nm))
  for (j in 1:4) {                       # same algebra as R/corner.R:48,55-56
    a <- cor[[j]][1]; ap <- cor[[j]][2]
    muAM  <- mu_oracle(a, M, C, s, tau, binY)
    ratio <- (qa(ap) / qa(a)) * (pa(a) / pa(ap))
    eta   <- eta_oracle(a, ap, C, s, tau, binY)
    phi[, j] <- (A == a) / pa(a) * ratio * (Y - muAM) +
                (A == ap) / pa(ap) * (muAM - eta) + eta
  }
  phi
}

## ---- exact truth (closed form; weak-ID grid's construction, with tau) --------
theta_exact <- function(a, ap, s, tau, binY) {
  b  <- 0.7 + tau * s * a
  mu <- 0.5 * s * a + 0.6 * s * b * ap
  if (!binY) return(mu)
  sdv <- sqrt((0.4 * b + 0.3)^2 + b^2)
  stats::integrate(function(x) plogis(x) * dnorm(x, mu, sdv),
                   mu - 12 * sdv, mu + 12 * sdv, rel.tol = 1e-10)$value
}
truth <- function(s, tau, binY) {
  th <- c(`11` = theta_exact(1, 1, s, tau, binY), `10` = theta_exact(1, 0, s, tau, binY),
          `01` = theta_exact(0, 1, s, tau, binY), `00` = theta_exact(0, 0, s, tau, binY))
  OE <- th["11"] - th["00"]; IDE <- th["10"] - th["00"]; IIE <- th["01"] - th["00"]
  list(W = unname((OE - IDE - IIE) / OE), OE = unname(OE), theta = th)
}

## ---- W and its IF SE from any influence matrix (R/gauge-pmed.R:262-272) -----
gauge_stats <- function(phi, zc = qnorm(0.975)) {
  n <- nrow(phi); t <- colMeans(phi)
  OE <- t["11"] - t["00"]; IDE <- t["10"] - t["00"]; IIE <- t["01"] - t["00"]
  W <- (OE - IDE - IIE) / OE
  pOE <- phi[, "11"] - phi[, "00"]; pIDE <- phi[, "10"] - phi[, "00"]
  pIIE <- phi[, "01"] - phi[, "00"]; pR <- pOE - pIDE - pIIE
  seW <- sd((pR - W * pOE) / OE) / sqrt(n)
  c(W = unname(W), seW = unname(seW), OE = unname(OE),
    oe_snr = unname(abs(OE) / (sd(pOE) / sqrt(n))),
    lo = unname(W - zc * seW), hi = unname(W + zc * seW))
}

## ---- cells -------------------------------------------------------------------
cells <- rbind(
  cbind(expand.grid(n = c(800L, 3000L), tau = c(0, 0.8), binY = c(FALSE, TRUE), s = 1),
        regime = "manuscript"),
  cbind(expand.grid(n = 800L, tau = 0.4, binY = c(FALSE, TRUE), s = c(0.2, 0.5)),
        regime = c("near-null", "near-null", "intermediate", "intermediate")))
cells$id <- seq_len(nrow(cells))

## ---- oracle wiring check: E[phi_oracle] must equal theta_exact ---------------
oracle_check <- function(s, tau, binY, n = 2e5L) {
  set.seed(1); d <- gen(n, s, tau, binY)
  th_hat <- colMeans(phi_oracle(d, s, tau, binY))
  th <- truth(s, tau, binY)$theta
  err <- max(abs(th_hat - th))
  if (err > 0.01) stop(sprintf("oracle EIF mean off truth by %.4f (s=%g tau=%g binY=%s)",
                               err, s, tau, binY))
  err
}

## ---- one cell -----------------------------------------------------------------
run_cell <- function(ci, nrep = NREP, R_part = R_PART) {
  cl <- cells[ci, ]; n <- cl$n; s <- cl$s; tau <- cl$tau; binY <- cl$binY
  tr <- truth(s, tau, binY)
  W_cf <- se_cf <- matrix(NA_real_, nrep, R_part)
  rows <- vector("list", nrep); nfail <- 0L; t0 <- proc.time()[["elapsed"]]
  base <- 1000000L * ci
  for (r in seq_len(nrep)) {
    seed <- base + r
    set.seed(seed); d <- gen(n, s, tau, binY)
    so <- gauge_stats(phi_oracle(d, s, tau, binY))
    ok <- TRUE
    for (p in seq_len(R_part)) {                 # same data, R_part fold partitions
      set.seed(seed * 100L + p)
      ph <- tryCatch(corner_fit(d, K, binY, "C")$phi, error = function(e) NULL)
      if (is.null(ph)) { ok <- FALSE; break }
      g <- gauge_stats(ph); W_cf[r, p] <- g["W"]; se_cf[r, p] <- g["seW"]
      if (p == 1L) g1 <- g
    }
    if (!ok) { nfail <- nfail + 1L; W_cf[r, ] <- se_cf[r, ] <- NA; next }
    rows[[r]] <- data.frame(
      cell = ci, seed = seed, trW = tr$W,
      W_cf1 = g1["W"], seW_an = g1["seW"], oe_snr = g1["oe_snr"],
      covW_an = tr$W >= g1["lo"] && tr$W <= g1["hi"],
      W_or = so["W"], seW_or = so["seW"], oe_snr_or = so["oe_snr"],
      covW_or = tr$W >= so["lo"] && tr$W <= so["hi"])
  }
  rows <- do.call(rbind, rows); keep <- !is.na(W_cf[, 1])
  W_cf <- W_cf[keep, , drop = FALSE]; se_cf <- se_cf[keep, , drop = FALSE]

  ## ---- the decomposition (variances across datasets) ----
  V_cf   <- var(W_cf[, 1])                              # as shipped: one partition
  Fv <- apply(W_cf, 1, var); F <- mean(Fv)              # E_data[Var_partition]
  V_avg  <- V_cf - F                                    # law of total variance
  V_avg2 <- var(rowMeans(W_cf)) - F / R_part            # direct cross-check
  V_or   <- var(rows$W_or)
  S_an   <- mean(rows$seW_an^2); S_or <- mean(rows$seW_or^2)
  G <- V_cf - S_an; N <- V_avg - V_or; R <- V_or - S_an
  elapsed <- proc.time()[["elapsed"]] - t0

  summary <- data.frame(
    cell = ci, regime = cl$regime, n = n, tau = tau, binY = binY, s = s,
    nrep = nrow(rows), nfail = nfail, R_part = R_part, trW = tr$W,
    oe_snr_med = median(rows$oe_snr), pct_oe_regular = mean(rows$oe_snr >= 2),
    ## positive controls -- must reproduce the published grid in manuscript cells
    se_ratio_an = sqrt(S_an / V_cf),                    # expect 0.65-0.83
    covW_an = mean(rows$covW_an),                       # expect 0.86-0.91
    bias_cf = mean(W_cf[, 1]) - tr$W, bias_or = mean(rows$W_or) - tr$W,
    ## the decomposition, on the SD scale and as shares of the gap
    empSD = sqrt(V_cf), sd_or = sqrt(V_or), seW_an = sqrt(S_an), seW_or = sqrt(S_or),
    gap = G, F_fold = F, N_nuis = N, R_rem = R,
    F_fold_med = median(Fv), F_tail = F / median(Fv),   # >>1: partition noise is heavy-tailed
    share_F = F / G, share_N = N / G, share_R = R / G,
    R_formula = V_or - S_or, R_plugin = S_or - S_an,    # sub-split of R
    V_avg_ltv = V_avg, V_avg_direct = V_avg2,           # must agree
    ## diagnostics of the oracle arm and of SE instability under re-partition
    se_ratio_or = sqrt(S_or / V_or),                    # ~1 if the IF formula is right
    covW_or = mean(rows$covW_or),
    cv_seW_fold = mean(apply(se_cf, 1, sd) / rowMeans(se_cf)),
    ## robust scale, for the near-null cells where variances are outlier-driven
    mad_cf = mad(W_cf[, 1]), mad_or = mad(rows$W_or),
    seW_an_med = median(rows$seW_an),
    se_ratio_an_robust = median(rows$seW_an) / mad(W_cf[, 1]),
    elapsed_s = round(elapsed, 1), row.names = NULL)
  out <- list(summary = summary, rows = rows, W_cf = W_cf, se_cf = se_cf, cell = cl)
  saveRDS(out, file.path(OUTDIR, sprintf("phi_decomp_cell_%02d.rds", ci)))
  summary
}

## ---- collate --------------------------------------------------------------------
collate <- function() {
  fs <- list.files(OUTDIR, "^phi_decomp_cell_\\d+\\.rds$", full.names = TRUE)
  if (!length(fs)) stop("no cell files in ", OUTDIR)
  S <- do.call(rbind, lapply(fs, function(f) readRDS(f)$summary))
  S <- S[order(S$cell), ]
  write.csv(S, file.path(OUTDIR, "phi_decomp_summary.csv"), row.names = FALSE)
  show <- c("cell", "regime", "n", "tau", "binY", "s", "nrep", "nfail",
            "se_ratio_an", "covW_an", "empSD", "seW_an", "sd_or", "seW_or",
            "share_F", "share_N", "share_R", "F_tail", "se_ratio_or", "cv_seW_fold",
            "oe_snr_med", "se_ratio_an_robust")
  num <- vapply(S[show], is.numeric, logical(1))
  S2 <- S[show]; S2[num] <- lapply(S2[num], function(x) signif(x, 3))
  print(S2, row.names = FALSE)
  invisible(S)
}

## ---- main ------------------------------------------------------------------------
if (COLLATE) {
  collate()
} else {
  sel <- flag("cell")
  ids <- if (!is.null(sel)) as.integer(strsplit(sel, ",")[[1]])
         else if (SMOKE) c(1L, 5L, 9L)
         else if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID")))
           as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
         else cells$id
  ## oracle wiring check once per distinct DGP among the selected cells
  dg <- unique(cells[ids, c("s", "tau", "binY")])
  for (i in seq_len(nrow(dg)))
    message(sprintf("oracle check s=%g tau=%g binY=%s: max|E[phi]-theta| = %.2e",
                    dg$s[i], dg$tau[i], dg$binY[i],
                    oracle_check(dg$s[i], dg$tau[i], dg$binY[i])))
  for (ci in ids) {
    message(sprintf("cell %d/%d: n=%d tau=%g binY=%s s=%g (%s) nrep=%d R_part=%d",
                    ci, nrow(cells), cells$n[ci], cells$tau[ci], cells$binY[ci],
                    cells$s[ci], cells$regime[ci], NREP, R_PART))
    s <- run_cell(ci)
    message(sprintf("   se_ratio_an=%.3f covW_an=%.3f | shares F=%.2f N=%.2f R=%.2f | %.0fs",
                    s$se_ratio_an, s$covW_an, s$share_F, s$share_N, s$share_R, s$elapsed_s))
  }
  if (length(ids) > 1L) collate()
}
