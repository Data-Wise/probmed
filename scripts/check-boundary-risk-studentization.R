#!/usr/bin/env Rscript
# T2 follow-up: do the 3 SD misses of check-boundary-risk-derivation.R come from
# studentization (the rule uses |d_hat|/sqrt(v_hat) and cor(d_hat, v_hat) is .72-.77)?
# Same pre-stated SD gate: predicted/observed SD in [.85, 1.15]. No fresh datasets.
#
# Arm P (requested): resample archived (d_hat, v_hat) PAIRS with replacement and apply the
#   rule. NOTE: this cannot fail. The archived estimate is a deterministic function of
#   the archived pair, so a pair bootstrap reproduces the observed SD up to resampling
#   noise. It is reported because it was asked for, not as evidence for the mechanism.
# Arm S (informative control): same rows, but v_hat is drawn from a different row than
#   d_hat, which breaks the d_hat/v_hat dependence and keeps both marginals. If
#   studentization causes the misses, Arm S should reproduce them and Arm P should not.
#
# Usage: Rscript scripts/check-boundary-risk-studentization.R <replicates.csv>
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, file.exists(args[1]))
r <- read.csv(args[1])
lam <- c(select95_fieller = stats::qnorm(.975), select995_fieller = stats::qnorm(.9975))

# Row-wise version of the point rule (each row has its own v_hat).
point_rows <- function(d1, d2, v11, v12, v22, lambda) {
  weak1 <- abs(d1) / sqrt(v11) <= lambda
  weak2 <- abs(d2) / sqrt(v22) <= lambda
  p <- sqrt(abs(d1)) / (sqrt(abs(d1)) + sqrt(abs(d2)))
  den1 <- d1^2 + 4 * v11; den2 <- d2^2 + 4 * v22
  r1 <- v11 / den1; r2 <- v22 / den2
  cross <- v12 * d1 * d2 / (den1 * den2)
  shift <- (-(1 + 2 * p) * r1 + (3 - 2 * p) * r2) / 8 - (1 - 2 * p) * cross / 4
  out <- stats::plogis(stats::qlogis(p) - shift)
  out[weak1 & !weak2] <- 0
  out[!weak1 & weak2] <- 1
  out
}

set.seed(92) # arbitrary prediction seed, not a study seed
B <- 2e5
rows <- list()
worst <- 0
for (cell in c(2, 4, 5, 6, 7)) {
  for (meth in names(lam)) {
    o <- r[r$cell == cell & r$method == meth, ]
    n <- nrow(o)
    # reimplementation check against archived estimates
    est <- point_rows(o$d1, o$d2, o$v11, o$v12, o$v22, lam[[meth]])
    worst <- max(worst, abs(est - o$estimate))
    osd <- sd(o$estimate); osnap <- mean(o$selected > 0)
    i <- sample.int(n, B, replace = TRUE)
    j <- sample.int(n, B, replace = TRUE)
    fp <- point_rows(o$d1[i], o$d2[i], o$v11[i], o$v12[i], o$v22[i], lam[[meth]])
    fs <- point_rows(o$d1[i], o$d2[i], o$v11[j], o$v12[j], o$v22[j], lam[[meth]])
    snap <- function(a, b, c) {
      w1 <- abs(o$d1[a]) / sqrt(o$v11[b]) <= lam[[meth]]
      w2 <- abs(o$d2[a]) / sqrt(o$v22[b]) <= lam[[meth]]
      mean(xor(w1, w2))
    }
    rows[[length(rows) + 1]] <- data.frame(
      cell, meth, obs_sd = osd, pair_sd = sd(fp), shuf_sd = sd(fs),
      pair_ratio = sd(fp) / osd, shuf_ratio = sd(fs) / osd,
      snap_obs = osnap, snap_pair = snap(i, i), snap_shuf = snap(i, j),
      cor_d1_v11 = cor(o$d1, o$v11)
    )
  }
}
cat(sprintf("Reimplementation vs archived estimates, max abs diff over all rows: %.2e\n\n", worst))
stopifnot(worst < 1e-9)
res <- do.call(rbind, rows)
res$pair_ok <- res$pair_ratio >= .85 & res$pair_ratio <= 1.15
res$shuf_ok <- res$shuf_ratio >= .85 & res$shuf_ratio <= 1.15
print(format(res, digits = 3), row.names = FALSE)
miss <- res$cell == 4 & res$meth == "select95_fieller" |
  res$cell %in% c(5, 7) & res$meth == "select995_fieller"
cat(sprintf("\nOriginal 3 misses (cell 4 select95; cells 5, 7 select995):\n"))
print(format(res[miss, c("cell", "meth", "obs_sd", "pair_sd", "shuf_sd", "pair_ok", "shuf_ok")], digits = 3),
  row.names = FALSE)
cat(sprintf("\nArm P (pair bootstrap) SD gate: %d/%d pass; on the 3 original misses: %d/3.\n",
  sum(res$pair_ok), nrow(res), sum(res$pair_ok[miss])))
cat(sprintf("Arm S (shuffled v_hat)  SD gate: %d/%d pass; on the 3 original misses: %d/3.\n",
  sum(res$shuf_ok), nrow(res), sum(res$shuf_ok[miss])))
