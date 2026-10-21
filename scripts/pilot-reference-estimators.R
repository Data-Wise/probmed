#!/usr/bin/env Rscript
# G8 pilot: achievable max RMSE of reference estimators at the gate grid of SPEC-boundary-gate-set-2026-10-21.md,
# on the PILOT seed block (1.5e9 + cell id * 10000 + replicate). REPORT-ONLY. Pre-registered BEFORE any run:
#   Reference set: raw plug-in P-hat (wb_binary$raw) and C(kappa, tau = .5), kappa in {0, 1, 2, 3.84} (T3 SPEC section 3),
#   no new tuning. psi_j = ((d_j^2 - kappa v_jj) + sqrt((d_j^2 - kappa v_jj)^2 + (tau v_jj)^2)) / 2, m_j = psi_j^(1/4), P = m1 / (m1 + m2).
#   Gate check per estimator and per (path, n): max over the 9 grid points of RMSE against the author's tolerance
#     .15 (direct, 500), .12 (indirect, 500), .125 (direct, 1000), .09 (indirect, 1000), using scripts/decide-gate-rmse.R:
#     PASS if RMSE + 2 MCSE <= tol, FAIL if RMSE - 2 MCSE > tol, else CONFIRM (the max over points of the verdict, FAIL > CONFIRM > PASS).
#   Pilot verdict (grill Q8, pre-committed): EXISTENCE if some reference estimator PASSes all four (path, n) gates;
#     NEGATIVE if every reference estimator FAILs at least one gate (then no T6 for these estimators, tolerances unchanged);
#     INDETERMINATE otherwise (the author decides; no tolerance changes either way).
# Usage: Rscript scripts/pilot-reference-estimators.R <out_dir> [R_per_cell=2000] [pr55_dir]
#   <pr55_dir> is REQUIRED (argument 3 or env PROBMED_PR55_DIR): the root of a PR #55 checkout, read only (wb_binary and the archive come from there).
# Control: wb_binary reproduces the archived raw_wald estimates of PR #55 cells 5, 6, 7 (selection-study seeds) to 1e-8.
args <- commandArgs(TRUE)
stopifnot(length(args) >= 1L)
out <- args[1]; dir.create(out, showWarnings = FALSE, recursive = TRUE)
R <- if (length(args) >= 2L) as.integer(args[2]) else 2000L
pr55 <- if (length(args) >= 3L) args[3] else Sys.getenv("PROBMED_PR55_DIR", "")
if (!nzchar(pr55) || !dir.exists(pr55)) stop("Give the PR #55 checkout root as the third argument or set PROBMED_PR55_DIR", call. = FALSE)

# quadrature rule copied from PR #55 R/wasserstein-observed.R (lines 3-14); the candidates script calls it as probmed:::.wasserstein_rule
.wasserstein_rule <- function(order, normal = FALSE) {
  j <- seq_len(order - 1L)
  off <- if (normal) sqrt(j) else j / sqrt(4 * j^2 - 1)
  mat <- matrix(0, order, order)
  mat[cbind(j, j + 1L)] <- off
  mat[cbind(j + 1L, j)] <- off
  eig <- eigen(mat, symmetric = TRUE)
  list(x = if (normal) eig$values else (eig$values + 1) / 2, w = eig$vectors[1L, ]^2)
}
src <- readLines(file.path(pr55, "scripts", "wasserstein-boundary-candidates.R"))
eval(parse(text = sub("probmed:::.wasserstein_rule", ".wasserstein_rule", src, fixed = TRUE)))
source("scripts/decide-gate-rmse.R") # defines decide() and mcse_rmse()

gen <- function(n, alpha, direct, inter) {
  C <- stats::rnorm(n); A <- stats::rbinom(n, 1, stats::plogis(.3 * C)); M <- alpha * A + .5 * C + stats::rnorm(n)
  Y <- stats::rbinom(n, 1, stats::plogis(-.5 + direct * A + .8 * M + inter * A * M + .3 * C))
  data.frame(A, M, Y, C)
}
corner_mean <- function(direct, alpha, inter, a, ap) {
  b <- .8 + inter * a; mu <- -.5 + direct * a + b * alpha * ap; sdv <- sqrt((.5 * b + .3)^2 + b^2)
  stats::integrate(function(z) stats::plogis(mu + sdv * z) * stats::dnorm(z), -Inf, Inf, rel.tol = 1e-11)$value
}
truth_p <- function(direct, alpha, inter) {
  d1 <- sqrt(abs(corner_mean(direct, alpha, inter, 1, 1) - corner_mean(direct, alpha, inter, 1, 0)))
  d2 <- sqrt(abs(corner_mean(direct, alpha, inter, 1, 0) - corner_mean(direct, alpha, inter, 0, 0)))
  d1 / (d1 + d2)
}
cfun <- function(d, v, kappa, tau = .5) {
  psi <- ((d^2 - kappa * diag(v)) + sqrt((d^2 - kappa * diag(v))^2 + (tau * diag(v))^2)) / 2
  m <- psi^.25; m[1] / sum(m)
}
ests <- c("plugin", "C0", "C1", "C2", "C3.84"); kap <- c(NA, 0, 1, 2, 3.84)

# control: reproduce archived raw_wald for cells 5-7, replicates 1-5
arch <- utils::read.csv(file.path(pr55, "docs/reports/data/wasserstein-selection-2026-10-20/replicates.csv"))
arch <- arch[arch$method == "raw_wald" & arch$cell %in% 5:7 & arch$replicate <= 5, ]
cs <- data.frame(cell = 5:7, alpha = c(.7, 0, .7), direct = c(0, .6, .05), inter = c(0, .2, 0))
dif <- vapply(seq_len(nrow(arch)), function(i) {
  s <- cs[cs$cell == arch$cell[i], ]; set.seed(arch$seed[i])
  abs(wb_binary(gen(500, s$alpha, s$direct, s$inter))$raw - arch$estimate[i])
}, 0)
cat(sprintf("Control: wb_binary reproduces %d archived raw_wald estimates, max abs difference %.2e (need < 1e-8)\n", length(dif), max(dif)))
stopifnot(max(dif) < 1e-8)
if (R < 100L) stop("R too small for a pilot run")

# gate grid (SPEC section 1) and cell ids 3-38
grid <- list(
  list(n = 500L, path = "direct", p = c(0, .00875, .0175, .02635, .0352, .05325, .0713, .1095, .1477)),
  list(n = 500L, path = "indirect", p = c(0, .0040, .0080, .0120, .0160, .0242, .0324, .05005, .0677)),
  list(n = 1000L, path = "direct", p = c(0, .00625, .0125, .01875, .025, .03775, .0505, .0777, .1049)),
  list(n = 1000L, path = "indirect", p = c(0, .00285, .0057, .00855, .0114, .01725, .0231, .03525, .0474)))
cells <- do.call(rbind, lapply(seq_along(grid), function(g) data.frame(n = grid[[g]]$n, path = grid[[g]]$path,
  param = grid[[g]]$p, id = 2L + 9L * (g - 1L) + seq_along(grid[[g]]$p))))
cells$cellid <- cells$id
run <- function(i) {
  s <- cells[i, ]
  alpha <- if (s$path == "direct") .7 else s$param; direct <- if (s$path == "direct") s$param else .6
  inter <- if (s$path == "direct") 0 else .2
  tp <- truth_p(direct, alpha, inter)
  err <- matrix(NA_real_, R, length(ests)); fail <- 0L
  for (r in seq_len(R)) {
    set.seed(1500000000L + s$cellid * 10000L + r)
    f <- tryCatch(wb_binary(gen(s$n, alpha, direct, inter)), error = function(e) NULL)
    if (is.null(f)) { fail <- fail + 1L; next }
    err[r, ] <- c(f$raw, vapply(kap[-1], function(k) tryCatch(cfun(f$d, f$v, k), error = function(e) NA_real_), 0)) - tp
  }
  ok <- stats::complete.cases(err)
  res <- data.frame(s[c("n", "path", "param", "id")], truth = tp, failed = fail + sum(!ok) - fail, estimator = ests,
    bias = colMeans(err[ok, , drop = FALSE]), rmse = sqrt(colMeans(err[ok, , drop = FALSE]^2)),
    mcse = apply(err[ok, , drop = FALSE], 2, mcse_rmse), R_ok = sum(ok), row.names = NULL)
  utils::write.csv(res, file.path(out, sprintf("pilot-cell-%02d.csv", s$id)), row.names = FALSE)
  message("cell ", s$id, " done")
  res
}
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(cells)), run, mc.cores = 4L))
utils::write.csv(res, file.path(out, "pilot-rmse.csv"), row.names = FALSE)

tol <- c(direct500 = .15, indirect500 = .12, direct1000 = .125, indirect1000 = .09)
rank <- c(PASS = 1L, CONFIRM = 2L, FAIL = 3L); nm <- names(rank)
verdict <- do.call(rbind, lapply(ests, function(e) {
  do.call(rbind, lapply(names(tol), function(k) {
    d <- res[res$estimator == e & res$path == sub("[0-9]+", "", k) & res$n == as.integer(sub("[a-z]+", "", k)), ]
    v <- mapply(decide, d$rmse, d$mcse, tol[[k]])
    data.frame(estimator = e, gate = k, tol = tol[[k]], max_rmse = max(d$rmse), at_param = d$param[which.max(d$rmse)],
      verdict = nm[max(rank[v])])
  }))
}))
cat("\nPilot gate verdicts (max over the 9 grid points; PASS < CONFIRM < FAIL):\n"); print(verdict, row.names = FALSE, digits = 3)
pass_all <- tapply(verdict$verdict == "PASS", verdict$estimator, all)
fail_any <- tapply(verdict$verdict == "FAIL", verdict$estimator, any)
cat(sprintf("\nPILOT VERDICT: %s\n", if (any(pass_all)) paste("EXISTENCE (passes all four gates):", paste(names(pass_all)[pass_all], collapse = ", "))
  else if (all(fail_any)) "NEGATIVE (every reference estimator FAILs at least one gate; no T6 for these estimators)" else "INDETERMINATE (author decides)"))
utils::write.csv(verdict, file.path(out, "pilot-verdicts.csv"), row.names = FALSE)
