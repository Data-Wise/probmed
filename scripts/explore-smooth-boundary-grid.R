#!/usr/bin/env Rscript
# T3 design exploration (pre-freeze, NOT a validation): predicted bias of the smooth
# candidate  m_j = psi(d_j; s_j)^(1/4),  psi = ((d^2 - kappa s^2) + sqrt((d^2 - kappa s^2)^2 + (tau s^2)^2)) / 2,
# P = m1 / (m1 + m2), under the Gaussian experiment d_hat ~ N(true d, Vbar) per cell.
# The Gaussian experiment predicted BIAS well (T2: 20/20) and SD poorly in rare-snap cells,
# so only bias and the zero-contrast bias coefficient are used to choose among designs.
# No fresh datasets; arbitrary prediction seed 96.
#
# Usage: Rscript scripts/explore-smooth-boundary-grid.R <replicates.csv>
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, file.exists(args[1]))
r <- read.csv(args[1])

cand <- function(d1, d2, v11, v22, kappa, tau) {
  psi <- function(d, v) {
    a <- d^2 - kappa * v
    (a + sqrt(a^2 + (tau * v)^2)) / 2
  }
  m1 <- psi(d1, v11)^(1 / 4); m2 <- psi(d2, v22)^(1 / 4)
  m1 / (m1 + m2)
}
# zero-contrast coefficient: E[m]/sqrt(s) for delta = 0 (the n^(-1/4) bias coefficient)
coef0 <- function(kappa, tau) {
  stats::integrate(function(z) {
    a <- z^2 - kappa
    ((a + sqrt(a^2 + tau^2)) / 2)^(1 / 4) * stats::dnorm(z)
  }, -Inf, Inf)$value
}
raw_coef <- 2^(1 / 4) * gamma(3 / 4) / sqrt(pi)

grid <- expand.grid(kappa = c(0, 1, 2, 3.84), tau = c(.5, 1, 2))
set.seed(96)
B <- 4e5
z0 <- matrix(stats::rnorm(B), ncol = 2); zz <- rbind(z0, -z0)
cells <- c(2, 4, 5, 6, 7)
out <- list()
for (cell in cells) {
  b <- r[r$cell == cell & r$method == "select95_fieller", ]
  d0 <- c(b$true_d1[1], b$true_d2[1]); truth <- b$truth[1]
  v <- c(mean(b$v11), mean(b$v12), mean(b$v22))
  x <- sweep(zz %*% chol(matrix(c(v[1], v[2], v[2], v[3]), 2)), 2, d0, "+")
  for (g in seq_len(nrow(grid))) {
    p <- cand(x[, 1], x[, 2], v[1], v[3], grid$kappa[g], grid$tau[g])
    out[[length(out) + 1]] <- data.frame(
      cell, kappa = grid$kappa[g], tau = grid$tau[g],
      pred_bias = mean(p) - truth, pred_sd = sd(p), bias_over_sd = (mean(p) - truth) / sd(p)
    )
  }
}
res <- do.call(rbind, out)
cat(sprintf("Raw plug-in zero-contrast coefficient E|Z|^(1/2) = %.3f\n", raw_coef))
cg <- cbind(grid, coef0 = mapply(coef0, grid$kappa, grid$tau))
cat("Candidate zero-contrast coefficient E[m]/sqrt(s) (smaller = less cusp bias, same n^(-1/4) rate):\n")
print(format(cg, digits = 3), row.names = FALSE)
w <- reshape(res[, c("cell", "kappa", "tau", "pred_bias")], idvar = c("kappa", "tau"),
  timevar = "cell", direction = "wide")
names(w) <- sub("pred_bias.", "bias_c", names(w))
cat("\nPredicted bias by cell (2,4 regular; 5 zero direct; 6 zero indirect; 7 weak direct):\n")
print(format(w, digits = 3), row.names = FALSE)
w2 <- reshape(res[, c("cell", "kappa", "tau", "bias_over_sd")], idvar = c("kappa", "tau"),
  timevar = "cell", direction = "wide")
names(w2) <- sub("bias_over_sd.", "b/sd_c", names(w2))
cat("\nPredicted |bias|/SD (frozen gate < .10; SD from the Gaussian experiment, treat as indicative):\n")
print(format(w2, digits = 3), row.names = FALSE)
