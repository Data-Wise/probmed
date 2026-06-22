#!/usr/bin/env Rscript
# Guard: the HPC-portable standalone (gauge_coverage_standalone.R) must stay
# numerically IDENTICAL to the packaged ward_residual() at the same seed.
# Run from the package root with devtools available:
#   Rscript inst/sim/check_standalone_equiv.R
# Exits non-zero if any difference exceeds tolerance.
suppressMessages(devtools::load_all(quiet = TRUE))
source("inst/sim/gauge_coverage_standalone.R", local = TRUE)

set.seed(42)
n <- 1000
C <- rnorm(n); A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
M <- 0.6 * A + 0.4 * C + rnorm(n)
Y <- 0.5 * A + 0.7 * M + 0.8 * A * M + 0.3 * C + rnorm(n)
d <- data.frame(A, M, Y, C)

tol <- 1e-8; bad <- FALSE
for (sm in c("analytic", "bootstrap")) {
  B <- if (sm == "bootstrap") 200L else 1L
  rp <- ward_residual(d, se_method = sm, B = B, fieller = FALSE, seed = 7L)
  rs <- gauge_core(d, se_method = sm, B = B, seed = 7L)
  dmax <- max(abs(rp@W - rs$W), abs(rp@p_med - rs$P),
              max(abs(rp@W_ci - rs$W_ci)), max(abs(rp@p_med_ci - rs$p_med_ci)))
  cat(sprintf("%-9s  max abs diff = %.2e  %s\n", sm, dmax,
              if (dmax <= tol) "OK" else "MISMATCH"))
  if (dmax > tol) bad <- TRUE
}
if (bad) stop("standalone diverged from ward_residual()") else cat("standalone == package\n")
