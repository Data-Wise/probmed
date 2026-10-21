#!/usr/bin/env Rscript
# T2: check the boundary point-risk derivation (docs/specs/DERIVE-boundary-point-risk-2026-10-21.md)
# against the archived selection study. No fresh datasets: the only RNG is a
# Gaussian-experiment prediction below (seed is arbitrary, not a study seed).
#
# Usage: Rscript scripts/check-boundary-risk-derivation.R <replicates.csv>
# <replicates.csv> is docs/reports/data/wasserstein-selection-2026-10-20/replicates.csv
# on branch codex/dev-main-review (PR #55); it is not on dev.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, file.exists(args[1]))
r <- read.csv(args[1])

# Independent reimplementation of the point rules (read from
# scripts/wasserstein-boundary-selection.R, then checked against archived rows).
point <- function(x, v, lambda) {
  weak <- sweep(abs(x), 2, sqrt(diag(v)), "/") <= lambda
  p <- sqrt(abs(x[, 1])) / rowSums(sqrt(abs(x)))
  den <- sweep(x^2, 2, 4 * diag(v), "+")
  rr <- sweep(1 / den, 2, diag(v), "*")
  cross <- v[1, 2] * x[, 1] * x[, 2] / (den[, 1] * den[, 2])
  shift <- (-(1 + 2 * p) * rr[, 1] + (3 - 2 * p) * rr[, 2]) / 8 -
    (1 - 2 * p) * cross / 4
  out <- stats::plogis(stats::qlogis(p) - shift)
  out[weak[, 1] & !weak[, 2]] <- 0
  out[!weak[, 1] & weak[, 2]] <- 1
  out
}
raw_p <- function(x) sqrt(abs(x[, 1])) / rowSums(sqrt(abs(x)))
lam <- c(select95_fieller = stats::qnorm(.975), select995_fieller = stats::qnorm(.9975),
         k4_fieller = -1)

# Step 0: reimplementation reproduces archived estimates row by row.
chk <- r[r$method %in% names(lam), ]
worst <- 0
for (m in names(lam)) {
  s <- chk[chk$method == m, ][1:300, ]
  for (i in seq_len(nrow(s))) {
    v <- matrix(c(s$v11[i], s$v12[i], s$v12[i], s$v22[i]), 2)
    worst <- max(worst, abs(point(matrix(c(s$d1[i], s$d2[i]), 1), v, lam[[m]]) - s$estimate[i]))
  }
}
cat(sprintf("Step 0: max |reimplemented - archived| over 900 rows = %.2e\n", worst))
stopifnot(worst < 1e-9)

# Step 1: Gaussian-experiment prediction. d_hat ~ N(true d, Vbar) per cell.
set.seed(91)
m <- 4e5
zz <- matrix(stats::rnorm(2 * m / 2), ncol = 2)
z <- rbind(zz, -zz) # antithetic
tol_bias <- function(mcse) max(4 * mcse, 0.03)
rows <- list()
for (cell in c(2, 4, 5, 6, 7)) {
  b <- r[r$cell == cell & r$method == "select95_fieller", ]
  d0 <- c(b$true_d1[1], b$true_d2[1]); truth <- b$truth[1]
  vb <- matrix(c(mean(b$v11), mean(b$v12), mean(b$v12), mean(b$v22)), 2)
  x <- sweep(z %*% chol(vb), 2, d0, "+")
  pred <- list(
    k4_fieller = point(x, vb, lam[["k4_fieller"]]),
    select95_fieller = point(x, vb, lam[["select95_fieller"]]),
    select995_fieller = point(x, vb, lam[["select995_fieller"]]),
    raw_wald = raw_p(x)
  )
  snap <- function(l) {
    w <- sweep(abs(x), 2, sqrt(diag(vb)), "/") <= l
    mean(xor(w[, 1], w[, 2]))
  }
  for (meth in names(pred)) {
    o <- r[r$cell == cell & r$method == meth, ]
    ob <- mean(o$estimate - truth); osd <- sd(o$estimate)
    omcse <- osd / sqrt(nrow(o))
    pb <- mean(pred[[meth]] - truth); psd <- sd(pred[[meth]])
    sr_obs <- if (meth %in% c("select95_fieller", "select995_fieller")) mean(o$selected > 0) else NA
    sr_pred <- if (is.na(sr_obs)) NA else snap(lam[[meth]])
    rows[[length(rows) + 1]] <- data.frame(
      cell, meth, obs_bias = ob, pred_bias = pb, obs_sd = osd, pred_sd = psd,
      snap_obs = sr_obs, snap_pred = sr_pred,
      bias_ok = abs(ob - pb) <= tol_bias(omcse),
      sd_ok = psd / osd >= .85 & psd / osd <= 1.15,
      snap_ok = is.na(sr_obs) | abs(sr_obs - sr_pred) <= .03
    )
  }
}
res <- do.call(rbind, rows)
print(format(res, digits = 3), row.names = FALSE)
cat(sprintf("\nGates: bias within max(4*MCSE, .03): %d/%d; SD ratio in [.85,1.15]: %d/%d; snap rate within .03: %d/%d\n",
  sum(res$bias_ok), nrow(res), sum(res$sd_ok), nrow(res),
  sum(res$snap_ok), nrow(res)))

# Step 2: positive controls. A wrong derivation must fail the same gates.
ctl_a <- res[res$meth == "raw_wald", ]
# (A) first-order delta method: regular plug-in has zero bias at every cell.
a_fail <- abs(ctl_a$obs_bias - 0) > pmax(4 * ctl_a$obs_sd / sqrt(2000), .03)
cat(sprintf("Control A (delta-method bias = 0) rejected in cells: %s (must include the 3 stress cells 5,6,7)\n",
  paste(ctl_a$cell[a_fail], collapse = ",")))
# (B) wrong threshold in the snap-rate formula (lambda - 1).
b <- r[r$cell == 2 & r$method == "select95_fieller", ]
d0 <- c(b$true_d1[1], b$true_d2[1])
vb <- matrix(c(mean(b$v11), mean(b$v12), mean(b$v12), mean(b$v22)), 2)
x <- sweep(z %*% chol(vb), 2, d0, "+")
w <- sweep(abs(x), 2, sqrt(diag(vb)), "/") <= (lam[["select95_fieller"]] - 1)
b_pred <- mean(xor(w[, 1], w[, 2]))
b_obs <- mean(b$selected > 0)
cat(sprintf("Control B (lambda-1 at cell 2): predicted snap %.4f vs observed %.4f -> %s\n",
  b_pred, b_obs, if (abs(b_pred - b_obs) > .03) "rejected (good)" else "NOT rejected (check is blind)"))

ok <- all(res$bias_ok, res$sd_ok, res$snap_ok) && all(c(5, 6, 7) %in% ctl_a$cell[a_fail]) &&
  abs(b_pred - b_obs) > .03
cat(sprintf("\nT2 overall: %s\n", if (ok) "PASS" else "FAIL"))
quit(status = if (ok) 0L else 1L)
