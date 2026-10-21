#!/usr/bin/env Rscript
# T2b: test the derivation itself (Codex Review A, finding 2). Tolerances below were fixed
# before the first run; none is changed after seeing results. Failures are reported, not tuned.
#
#  (a) RATE (D1). Raw plug-in P_hat = sqrt|X1| / (sqrt|X1| + sqrt|X2|), X1 ~ N(d1, s1^2),
#      X2 ~ N(0, s2^2) independent, d1 = .1123 (cell 5), s1/s2 = .522 (cell 5), s2 =
#      .0428 / 2^k, k = 0..8 (halving s = quadrupling n). Bias computed by deterministic
#      2D midpoint quadrature (no RNG).
#      A1: log-log slope of |bias| on s2 over the 4 smallest s2 lies in [.45, .55].
#      A2: bias / (-E|Z|^(1/2) sqrt(s2) / sqrt(d1)) at the smallest s2 lies in [.93, 1.07].
#      A3 (planted control): the same slope for a REGULAR configuration (d2 = .1191, same s
#          grid) must lie OUTSIDE [.45, .55], i.e. the check can tell n^-1/4 from smooth.
#  (b) SNAP FORMULA (D2). For the 5 binary cells and lambda in {qnorm(.975), qnorm(.9975)}:
#      B1: independence-product snap probability within .01 of the exact correlated value
#          (if not, D2 must be restated with the correlation);
#      B2: exact correlated value within 4 Monte Carlo SE (+1e-4) of a 4e6-draw simulation
#          of the same Gaussian experiment;
#      B3: exact correlated value within .03 of the archived snap rate (as in T2).
#  (c) RISK DECOMPOSITION. RMSE^2 = sum_b pi_b [(m_b - truth)^2 + v_b] over branches
#      snap-to-0, snap-to-1, no-snap (k4). pi_b exact from (b); no-snap moments from the
#      Gaussian experiment. C1: predicted RMSE / observed RMSE in [.85, 1.15] for the 10
#      cell-method pairs. This is the Gaussian experiment organized by branch: it tests the
#      decomposition and the closed-form branch probabilities, not the Gaussian assumption.
#
# Usage: Rscript scripts/check-boundary-risk-rate.R <replicates.csv>
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, file.exists(args[1]))
r <- read.csv(args[1])
ok <- list()

# ---- (a) rate -------------------------------------------------------------
bias_raw <- function(d1, d2, s1, s2, truth, h = .006, lim = 9) {
  t <- seq(-lim + h / 2, lim - h / 2, by = h)
  w <- stats::dnorm(t) * h
  x1 <- sqrt(abs(d1 + s1 * t)); x2 <- sqrt(abs(d2 + s2 * t))
  total <- 0
  for (i in seq_along(t)) total <- total + w[i] * sum(w * (x1[i] / (x1[i] + x2)))
  total / (sum(w)^2) - truth
}
d1 <- .1123045; s2g <- .0428 / 2^(0:8); ratio <- .522
bz <- vapply(s2g, function(s) bias_raw(d1, 0, ratio * s, s, 1), 0)
sl <- function(b, s) unname(coef(lm(log(abs(b)) ~ log(s)))[2])
last4 <- 6:9
slope_zero <- sl(bz[last4], s2g[last4])
lead <- -(2^(1 / 4) * gamma(3 / 4) / sqrt(pi)) * sqrt(s2g[9]) / sqrt(d1)
coef_ratio <- bz[9] / lead
breg <- vapply(s2g, function(s) bias_raw(.1287156, .1190775, ratio * s, s, .5097276), 0)
slope_reg <- sl(breg[last4], s2g[last4])
cat("(a) RATE\n")
print(data.frame(s2 = s2g, bias_zero_contrast = bz, bias_regular = breg), digits = 4)
cat(sprintf("A1 slope (zero contrast, 4 smallest s) = %.3f, gate [.45,.55] -> %s\n",
  slope_zero, if (slope_zero >= .45 && slope_zero <= .55) "PASS" else "FAIL"))
cat(sprintf("A2 bias/leading-term at smallest s = %.3f, gate [.93,1.07] -> %s\n",
  coef_ratio, if (coef_ratio >= .93 && coef_ratio <= 1.07) "PASS" else "FAIL"))
cat(sprintf("A3 control: regular-configuration slope = %.3f, must be outside [.45,.55] -> %s\n\n",
  slope_reg, if (slope_reg < .45 || slope_reg > .55) "PASS (check discriminates)" else "FAIL (blind)"))
ok$A1 <- slope_zero >= .45 && slope_zero <= .55
ok$A2 <- coef_ratio >= .93 && coef_ratio <= 1.07
ok$A3 <- slope_reg < .45 || slope_reg > .55

# ---- (b) and (c) per cell ---------------------------------------------------
point_vec <- function(x1, x2, v11, v12, v22, lambda) {
  weak1 <- abs(x1) / sqrt(v11) <= lambda; weak2 <- abs(x2) / sqrt(v22) <= lambda
  p <- sqrt(abs(x1)) / (sqrt(abs(x1)) + sqrt(abs(x2)))
  den1 <- x1^2 + 4 * v11; den2 <- x2^2 + 4 * v22
  r1 <- v11 / den1; r2 <- v22 / den2
  cross <- v12 * x1 * x2 / (den1 * den2)
  shift <- (-(1 + 2 * p) * r1 + (3 - 2 * p) * r2) / 8 - (1 - 2 * p) * cross / 4
  list(k4 = stats::plogis(stats::qlogis(p) - shift), weak1 = weak1, weak2 = weak2)
}
set.seed(95) # arbitrary prediction seed, not a study seed
B <- 4e6
z0 <- matrix(stats::rnorm(B), ncol = 2); zz <- rbind(z0, -z0)
lams <- c(select95_fieller = stats::qnorm(.975), select995_fieller = stats::qnorm(.9975))
rows <- list()
for (cell in c(2, 4, 5, 6, 7)) {
  b <- r[r$cell == cell & r$method == "select95_fieller", ]
  d0 <- c(b$true_d1[1], b$true_d2[1]); truth <- b$truth[1]
  v11 <- mean(b$v11); v12 <- mean(b$v12); v22 <- mean(b$v22)
  s1 <- sqrt(v11); s2 <- sqrt(v22); rho <- v12 / (s1 * s2)
  mu <- d0 / c(s1, s2)
  x <- sweep(zz %*% chol(matrix(c(v11, v12, v12, v22), 2)), 2, d0, "+")
  pt <- point_vec(x[, 1], x[, 2], v11, v12, v22, -1)
  for (meth in names(lams)) {
    lam <- lams[[meth]]
    pw <- function(j) stats::pnorm(lam - mu[j]) - stats::pnorm(-lam - mu[j])
    lo <- c(-lam - mu[1], -lam - mu[2]); hi <- c(lam - mu[1], lam - mu[2])
    inner <- function(u) {
      stats::dnorm(u) * (stats::pnorm((hi[2] - rho * u) / sqrt(1 - rho^2)) -
        stats::pnorm((lo[2] - rho * u) / sqrt(1 - rho^2)))
    }
    p11 <- stats::integrate(inner, lo[1], hi[1], rel.tol = 1e-10)$value
    pi0 <- pw(1) - p11; pi1 <- pw(2) - p11; pin <- 1 - pi0 - pi1
    snap_exact <- pi0 + pi1
    snap_indep <- pw(1) + pw(2) - 2 * pw(1) * pw(2)
    w1 <- abs(x[, 1]) / s1 <= lam; w2 <- abs(x[, 2]) / s2 <= lam
    snap_mc <- mean(xor(w1, w2))
    mcse <- sqrt(snap_exact * (1 - snap_exact) / nrow(x))
    ns <- !(w1 & !w2) & !(!w1 & w2)
    mn <- mean(pt$k4[ns]); vn <- mean((pt$k4[ns] - mn)^2)
    pred_mse <- pi0 * (0 - truth)^2 + pi1 * (1 - truth)^2 + pin * ((mn - truth)^2 + vn)
    o <- r[r$cell == cell & r$method == meth, ]
    obs_rmse <- sqrt(mean((o$estimate - truth)^2))
    rows[[length(rows) + 1]] <- data.frame(
      cell, meth, snap_exact, snap_indep, snap_mc, snap_obs = mean(o$selected > 0),
      B1 = abs(snap_indep - snap_exact) <= .01,
      B2 = abs(snap_exact - snap_mc) <= 4 * mcse + 1e-4,
      B3 = abs(snap_exact - mean(o$selected > 0)) <= .03,
      pred_rmse = sqrt(pred_mse), obs_rmse, C1 = sqrt(pred_mse) / obs_rmse >= .85 &
        sqrt(pred_mse) / obs_rmse <= 1.15
    )
  }
}
res <- do.call(rbind, rows)
cat("(b) SNAP FORMULA and (c) RISK DECOMPOSITION\n")
print(format(res, digits = 3), row.names = FALSE)
for (g in c("B1", "B2", "B3", "C1")) cat(sprintf("%s: %d/%d pass\n", g, sum(res[[g]]), nrow(res)))
failed <- c(names(ok)[!unlist(ok)], c("B1", "B2", "B3", "C1")[
  vapply(c("B1", "B2", "B3", "C1"), function(g) !all(res[[g]]), NA)])
cat(sprintf("\nT2b overall: %s%s\n", if (length(failed) == 0) "PASS" else "FAIL",
  if (length(failed)) paste0(" (failed: ", paste(failed, collapse = ", "), ")") else ""))
quit(status = if (length(failed) == 0) 0L else 1L)
