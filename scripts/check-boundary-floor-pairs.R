#!/usr/bin/env Rscript
# G2b: two-point floors for ALL node pairs (not only boundary-to-node) on each ladder.
# Usage: Rscript scripts/check-boundary-floor-pairs.R [n]   (n in 500, 1000; default 500)
# Grid per path = boundary (0), the four TV-matched ladder positions from DERIVE-boundary-floor
# section 2/6, the cell-7-like point .05 (direct) / .05 (indirect), and midpoints of consecutive nodes.
# For each pair (i < j): data from law i, exact likelihood ratio of law j to law i,
#   TV_ij = E_i[(1 - Q_j/P_i)_+], gap_ij = |P(theta_j) - P(theta_i)|, floor_ij = gap_ij sqrt((1 - TV_ij) / 8)
# (the Le Cam bound of check-boundary-floor-ladder.R, now for any pair).
# Controls (script stops if any fails): (i) identical-law TV is 0; (ii) the boundary pair at the shared
# point .05 reproduces the boundary-to-node TV of floor-ladder-n<n>.csv within 4 combined MC SE.
args <- commandArgs(TRUE)
n <- if (length(args)) as.integer(args[1]) else 500L
stopifnot(n %in% c(500L, 1000L))
set.seed(20261023 + n)
R <- 10000L
pos <- list(
  `500`  = list(direct = c(.0175, .0352, .0713, .1477), indirect = c(.0080, .0160, .0324, .0677)),
  `1000` = list(direct = c(.0125, .0250, .0505, .1049), indirect = c(.0057, .0114, .0231, .0474)))[[as.character(n)]]
mkgrid <- function(p, extra) {
  g <- sort(unique(c(0, p, extra)))
  sort(unique(c(g, (g[-1] + g[-length(g)]) / 2)))
}
gd <- mkgrid(pos$direct, .05); ga <- mkgrid(pos$indirect, .05)

corner_mean <- function(direct, alpha, inter, a, ap) {
  b <- .8 + inter * a
  mu <- -.5 + direct * a + b * alpha * ap
  sdv <- sqrt((.5 * b + .3)^2 + b^2)
  stats::integrate(function(z) stats::plogis(mu + sdv * z) * stats::dnorm(z), -Inf, Inf, rel.tol = 1e-11)$value
}
truth_p <- function(direct, alpha, inter) {
  c11 <- corner_mean(direct, alpha, inter, 1, 1); c10 <- corner_mean(direct, alpha, inter, 1, 0)
  c00 <- corner_mean(direct, alpha, inter, 0, 0)
  d1 <- sqrt(abs(c11 - c10)); d2 <- sqrt(abs(c10 - c00))
  d1 / (d1 + d2)
}
# one base law i: returns TV of every law j (vector over grid) from datasets drawn under law i
tv_from_base_direct <- function(i, grid) {
  v <- vapply(seq_len(R), function(r) {
    C <- stats::rnorm(n); A <- stats::rbinom(n, 1, stats::plogis(.3 * C))
    M <- .7 * A + .5 * C + stats::rnorm(n)
    e0 <- -.5 + .8 * M + .3 * C
    Y <- stats::rbinom(n, 1, stats::plogis(e0 + grid[i] * A))
    li <- stats::dbinom(Y, 1, stats::plogis(e0 + grid[i] * A), log = TRUE)
    vapply(grid, function(d) max(0, 1 - exp(sum(stats::dbinom(Y, 1, stats::plogis(e0 + d * A), log = TRUE) - li))), 0)
  }, numeric(length(grid)))
  cbind(tv = rowMeans(v), se = apply(v, 1, stats::sd) / sqrt(R))
}
tv_from_base_alpha <- function(i, grid) {
  v <- vapply(seq_len(R), function(r) {
    C <- stats::rnorm(n); A <- stats::rbinom(n, 1, stats::plogis(.3 * C))
    M <- grid[i] * A + .5 * C + stats::rnorm(n)
    li <- stats::dnorm(M, grid[i] * A + .5 * C, 1, log = TRUE)
    vapply(grid, function(al) max(0, 1 - exp(sum(stats::dnorm(M, al * A + .5 * C, 1, log = TRUE) - li))), 0)
  }, numeric(length(grid)))
  cbind(tv = rowMeans(v), se = apply(v, 1, stats::sd) / sqrt(R))
}
pairs_for <- function(label, grid, tvfun, truthfun) {
  th <- vapply(grid, truthfun, 0)
  rows <- list()
  for (i in seq_along(grid)) {
    m <- tvfun(i, grid)
    if (m[i, "tv"] != 0) stop("control failed: identical-law TV is not 0")
    for (j in seq_along(grid)) if (j > i) {
      gap <- abs(th[j] - th[i])
      rows[[length(rows) + 1L]] <- data.frame(path = label, a = grid[i], b = grid[j], gap = gap,
        tv = m[j, "tv"], se = m[j, "se"], floor_rmse = gap * sqrt((1 - m[j, "tv"]) / 8))
    }
  }
  do.call(rbind, rows)
}
pd <- pairs_for("direct", gd, tv_from_base_direct, function(d) truth_p(d, .7, 0))
pa <- pairs_for("indirect", ga, tv_from_base_alpha, function(al) truth_p(.6, al, .2))
all <- rbind(pd, pa)

ref <- utils::read.csv(sprintf("docs/specs/data/floor-ladder-n%d.csv", n))
ctl <- function(path, p) {
  r <- all[all$path == path & all$a == 0 & abs(all$b - p) < 1e-12, ]
  q <- ref[ref$path == path & abs(ref$param - p) < 1e-12, ]
  abs(r$tv - q$tv) < 4 * sqrt(r$se^2 + q$se^2)
}
okc <- ctl("direct", .05) && ctl("indirect", .05)
cat(sprintf("Controls: identical-law TV 0 (checked per base); boundary pair at .05 matches floor-ladder-n%d.csv: %s\n", n, okc))
stopifnot(okc)

cat(sprintf("\nn = %d, R = %d datasets per base law. Largest floors by pair type:\n", n, R))
for (path in c("direct", "indirect")) {
  d <- all[all$path == path, ]
  bd <- d[d$a == 0, ]; nd <- d[d$a > 0, ]
  cat(sprintf("  %-8s boundary-to-node max %.4f (b = %.4f);  node-to-node max %.4f (a = %.4f, b = %.4f)\n",
    path, max(bd$floor_rmse), bd$b[which.max(bd$floor_rmse)],
    max(nd$floor_rmse), nd$a[which.max(nd$floor_rmse)], nd$b[which.max(nd$floor_rmse)]))
}
utils::write.csv(all, sprintf("docs/specs/data/floor-pairs-n%d.csv", n), row.names = FALSE)
cat("\nAll pairs written to docs/specs/data/floor-pairs-n", n, ".csv\n", sep = "")
