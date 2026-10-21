#!/usr/bin/env Rscript
# T2 follow-up 2: cell 4 select95 SD miss (fixed-V Gaussian predicts .0905, observed .0770,
# ratio 1.175; shuffling v_hat did NOT reproduce it, ratio 1.04). Swap one ingredient at a
# time from the fixed-V Gaussian toward the archive. Same pre-stated SD gate:
# predicted/observed SD in [.85, 1.15]. No fresh datasets; arbitrary prediction seed 93.
#
#   G0 d_hat ~ N(true d, Vbar), v_hat = Vbar            (the original prediction)
#   G1 d_hat resampled from archive (joint rows), v_hat = Vbar   -> tests d_hat shape/scale
#   G2 d_hat ~ N(true d, Vbar), v_hat resampled from archive     -> tests v_hat variability
#   G3 d_hat ~ N(true d, Cov_emp(d_hat)), v_hat = Vbar           -> tests Vbar vs empirical cov
#   G4 d_hat resampled, v_hat resampled independently             (= Arm S of the earlier check)
#
# Usage: Rscript scripts/check-boundary-risk-cell4.R <replicates.csv>
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, file.exists(args[1]))
r <- read.csv(args[1])
lam <- stats::qnorm(.975)

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
skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurt <- function(x) mean((x - mean(x))^4) / sd(x)^4 - 3

o <- r[r$cell == 4 & r$method == "select95_fieller", ]
d0 <- c(o$true_d1[1], o$true_d2[1])
vb <- c(mean(o$v11), mean(o$v12), mean(o$v22))
Vb <- matrix(c(vb[1], vb[2], vb[2], vb[3]), 2)
ce <- cov(cbind(o$d1, o$d2))
cat("Cell 4 (binary n = 1000), archived d_hat vs fixed-V Gaussian assumptions\n")
cat(sprintf("  mean d_hat  : %.5f %.5f   true %.5f %.5f\n", mean(o$d1), mean(o$d2), d0[1], d0[2]))
cat(sprintf("  sd d_hat    : %.5f %.5f   sqrt(Vbar diag) %.5f %.5f  ratio %.3f %.3f\n",
  sd(o$d1), sd(o$d2), sqrt(vb[1]), sqrt(vb[3]), sd(o$d1) / sqrt(vb[1]), sd(o$d2) / sqrt(vb[3])))
cat(sprintf("  cov d_hat   : %.3e vs Vbar12 %.3e (corr %.3f vs %.3f)\n",
  ce[1, 2], vb[2], cor(o$d1, o$d2), vb[2] / sqrt(vb[1] * vb[3])))
cat(sprintf("  skew/ex.kurt: d1 %.3f/%.3f  d2 %.3f/%.3f  (Gaussian 0/0)\n",
  skew(o$d1), kurt(o$d1), skew(o$d2), kurt(o$d2)))
cat(sprintf("  v_hat cv    : v11 %.3f  v22 %.3f ; cor(v11,v22) %.3f\n",
  sd(o$v11) / mean(o$v11), sd(o$v22) / mean(o$v22), cor(o$v11, o$v22)))
z1 <- (o$d1 - d0[1]) / sqrt(o$v11); z2 <- (o$d2 - d0[2]) / sqrt(o$v22)
cat(sprintf("  studentized (d_hat - true)/se_hat: mean %.3f %.3f  sd %.3f %.3f\n\n",
  mean(z1), mean(z2), sd(z1), sd(z2)))

set.seed(93)
B <- 4e5
n <- nrow(o)
gauss <- function(S) {
  zz <- matrix(stats::rnorm(B), ncol = 2)
  sweep(rbind(zz, -zz) %*% chol(S), 2, d0, "+")
}
obs_sd <- sd(o$estimate); obs_snap <- mean(o$selected > 0)
gate <- function(x) x >= .85 & x <= 1.15
out <- list()
add <- function(label, est, snap) {
  out[[length(out) + 1]] <<- data.frame(
    variant = label, pred_sd = sd(est), ratio = sd(est) / obs_sd,
    snap = snap, pass = gate(sd(est) / obs_sd))
}
snap_of <- function(d1, d2, v11, v22) mean(xor(abs(d1) / sqrt(v11) <= lam, abs(d2) / sqrt(v22) <= lam))

x <- gauss(Vb)
add("G0 Gaussian d, fixed Vbar", point_rows(x[, 1], x[, 2], vb[1], vb[2], vb[3], lam),
  snap_of(x[, 1], x[, 2], vb[1], vb[3]))
i <- sample.int(n, B, replace = TRUE)
add("G1 archived d, fixed Vbar", point_rows(o$d1[i], o$d2[i], vb[1], vb[2], vb[3], lam),
  snap_of(o$d1[i], o$d2[i], vb[1], vb[3]))
j <- sample.int(n, B, replace = TRUE)
x <- gauss(Vb)
add("G2 Gaussian d, archived v_hat", point_rows(x[, 1], x[, 2], o$v11[j], o$v12[j], o$v22[j], lam),
  snap_of(x[, 1], x[, 2], o$v11[j], o$v22[j]))
x <- gauss(ce)
x[, 1] <- x[, 1] - 0 # centered at true d by gauss(); empirical covariance supplies the spread
add("G3 Gaussian d with empirical Cov(d_hat), fixed Vbar", point_rows(x[, 1], x[, 2], vb[1], vb[2], vb[3], lam),
  snap_of(x[, 1], x[, 2], vb[1], vb[3]))
add("G4 archived d, independent archived v_hat", point_rows(o$d1[i], o$d2[i], o$v11[j], o$v12[j], o$v22[j], lam),
  snap_of(o$d1[i], o$d2[i], o$v11[j], o$v22[j]))
res <- do.call(rbind, out)
cat(sprintf("Observed SD %.4f, observed snap rate %.4f\n", obs_sd, obs_snap))
print(format(res, digits = 3), row.names = FALSE)
cat(sprintf("\nVariants within the gate: %s\n", paste(res$variant[res$pass], collapse = "; ")))

# Sampling-noise check. The archive has only n = 2000 replicates and ~38 snapped fits, so its
# own SD is noisy. Draw 4000 studies of n = 2000 from G0 and report how often the study SD is
# at or below the observed one, and the 95% band of the study SD / population SD.
set.seed(94)
R <- 4000
sds <- numeric(R)
for (k in seq_len(R)) {
  zz <- matrix(stats::rnorm(2 * n), ncol = 2)
  xx <- sweep(zz %*% chol(Vb), 2, d0, "+")
  sds[k] <- sd(point_rows(xx[, 1], xx[, 2], vb[1], vb[2], vb[3], lam))
}
cat(sprintf("\nSampling noise (G0, %d studies of n = %d): study SD 2.5/50/97.5%% = %.4f / %.4f / %.4f\n",
  R, n, quantile(sds, .025), quantile(sds, .5), quantile(sds, .975)))
cat(sprintf("P(study SD <= observed %.4f) = %.4f ; study SD / G0 population SD in [.85,1.15]: %.3f of studies\n",
  obs_sd, mean(sds <= obs_sd), mean(sds / res$pred_sd[1] >= .85 & sds / res$pred_sd[1] <= 1.15)))
