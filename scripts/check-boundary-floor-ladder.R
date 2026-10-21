#!/usr/bin/env Rscript
# G2: two-point (Le Cam) floor and weak-signal ladder positions for the binary Wasserstein
# proportion mediated, n = 500. Needs no candidate estimator and no archived estimates.
#
# Estimand (as in scripts/validate-wasserstein-joint-contrast.R on the PR #55 branch):
#   P = sqrt|d1| / (sqrt|d1| + sqrt|d2|), d1 = c(1,1) - c(1,0), d2 = c(1,0) - c(0,0),
#   c(a, a') = E plogis(mu + sd Z), mu = -.5 + direct a + b alpha a', b = .8 + interaction a,
#   sd = sqrt((.5 b + .3)^2 + b^2).
# Two ladders from the boundary law, each a one-parameter path in the DGP:
#   "direct":   direct in grid, alpha = .7, interaction = 0 (cells 5, 7 family; boundary direct = 0, P = 1)
#   "indirect": alpha  in grid, direct = .6, interaction = .2 (cell 6 family; boundary alpha = 0, P = 0)
# Data-law TV to the boundary law at n = 500: TV = E_boundary[(1 - Q/P)_+] with Q/P the exact
# likelihood ratio, estimated by Monte Carlo with common datasets across the grid.
# Floor: for ANY estimator T and two data laws with truths theta0, theta1 (gap D) and TV tau,
#   max_i E_i (T - theta_i)^2 >= D^2 (1 - tau) / 8   (test of the nearer truth; each error costs >= (D/2)^2),
#   so max RMSE >= D sqrt((1 - tau) / 8). Reported per ladder point.
# Controls (must pass or the script stops): shift 0 gives TV = 0; the direct shift .05 reproduces the
# cell 5 / cell 7 TV .1404 of scripts/check-gate-feasibility-tv.R within 4 MC SE; a large shift gives TV > .9.
# Usage: Rscript scripts/check-boundary-floor-ladder.R [n]   (n defaults to 500; the .1404 control applies at 500 only)
args <- commandArgs(TRUE)
n <- if (length(args)) as.integer(args[1]) else 500L
stopifnot(n %in% c(500L, 1000L))
set.seed(if (n == 500L) 20261022 else 20261022 + n) # n = 500 keeps the seed of the archived output
R <- 20000L

corner_mean <- function(direct, alpha, inter, a, ap) {
  b <- .8 + inter * a
  mu <- -.5 + direct * a + b * alpha * ap
  sdv <- sqrt((.5 * b + .3)^2 + b^2)
  stats::integrate(function(z) stats::plogis(mu + sdv * z) * stats::dnorm(z),
    -Inf, Inf, rel.tol = 1e-11)$value
}
truth_p <- function(direct, alpha, inter) {
  c11 <- corner_mean(direct, alpha, inter, 1, 1)
  c10 <- corner_mean(direct, alpha, inter, 1, 0)
  c00 <- corner_mean(direct, alpha, inter, 0, 0)
  d1 <- sqrt(abs(c11 - c10)); d2 <- sqrt(abs(c10 - c00))
  if (d1 + d2 == 0) NA_real_ else d1 / (d1 + d2)
}

# direct ladder: Y | A, M, C logit shifts by `direct * A`; the M and C parts of the law are unchanged
tv_direct <- function(grid) {
  v <- vapply(seq_len(R), function(i) {
    C <- stats::rnorm(n); A <- stats::rbinom(n, 1, stats::plogis(.3 * C))
    M <- .7 * A + .5 * C + stats::rnorm(n)
    e0 <- -.5 + .8 * M + .3 * C
    Y <- stats::rbinom(n, 1, stats::plogis(e0))
    vapply(grid, function(d) {
      ll <- sum(stats::dbinom(Y, 1, stats::plogis(e0 + d * A), log = TRUE) -
        stats::dbinom(Y, 1, stats::plogis(e0), log = TRUE))
      max(0, 1 - exp(ll))
    }, 0)
  }, numeric(length(grid)))
  cbind(tv = rowMeans(v), se = apply(v, 1, stats::sd) / sqrt(R))
}
# indirect ladder: M | A, C mean shifts by `alpha * A`; Y | A, M, C is unchanged, so only the M part of the law moves
tv_alpha <- function(grid) {
  v <- vapply(seq_len(R), function(i) {
    C <- stats::rnorm(n); A <- stats::rbinom(n, 1, stats::plogis(.3 * C))
    M <- 0 * A + .5 * C + stats::rnorm(n)
    vapply(grid, function(al) {
      ll <- sum(stats::dnorm(M, al * A + .5 * C, 1, log = TRUE) - stats::dnorm(M, .5 * C, 1, log = TRUE))
      max(0, 1 - exp(ll))
    }, 0)
  }, numeric(length(grid)))
  cbind(tv = rowMeans(v), se = apply(v, 1, stats::sd) / sqrt(R))
}

report <- function(label, grid, tvm, truth_fun, p0) {
  p <- vapply(grid, truth_fun, 0)
  gap <- abs(p - p0)
  out <- data.frame(param = grid, theta = p, gap = gap, tv = tvm[, "tv"], se = tvm[, "se"],
    floor_rmse = gap * sqrt((1 - tvm[, "tv"]) / 8))
  cat(sprintf("\n%s ladder (n = %d, R = %d datasets per row; boundary truth %.4f)\n", label, n, R, p0))
  print(out, row.names = FALSE, digits = 4)
  out
}

dgrid <- c(0, .01, .02, .05, .10, .15, .25, .5, 1)
agrid <- c(0, .02, .05, .10, .15, .25, .5, 1)
td <- tv_direct(dgrid)
ta <- tv_alpha(agrid)

# controls
ok0 <- td[1, "tv"] == 0 && ta[1, "tv"] == 0
i05 <- which(dgrid == .05)
ok05 <- if (n == 500L) abs(td[i05, "tv"] - .1404) < 4 * sqrt(td[i05, "se"]^2 + .0012^2) else td[i05, "tv"] > .1418 # n = 1000: more data, TV must exceed the n = 500 value
okbig <- td[length(dgrid), "tv"] > .9 && ta[length(agrid), "tv"] > .9
cat(sprintf("Controls: shift 0 gives TV 0: %s; direct .05 gives %.4f (n = 500 value .1404/.1418; at n = 1000 must exceed it): %s; largest shifts give TV > .9: %s\n",
  ok0, td[i05, "tv"], ok05, okbig))
stopifnot(ok0, ok05, okbig)

rd <- report("direct", dgrid, td, function(d) truth_p(d, .7, 0), truth_p(0, .7, 0))
ra <- report("indirect", agrid, ta, function(al) truth_p(.6, al, .2), truth_p(.6, 0, .2))

# ladder positions at target data-law TV (linear interpolation of the param on TV; bracketed grid only)
targets <- c(.05, .10, .20, .40)
pos <- function(r) vapply(targets, function(t) {
  if (t > max(r$tv)) return(NA_real_)
  stats::approx(r$tv, r$param, xout = t, ties = "ordered")$y
}, 0)
csv <- rbind(cbind(path = "direct", rd), cbind(path = "indirect", ra))
utils::write.csv(csv, sprintf("docs/specs/data/floor-ladder-n%d.csv", n), row.names = FALSE)
cat("\nLadder positions at data-law TV to the boundary law:\n")
print(data.frame(target_tv = targets, direct = pos(rd), indirect = pos(ra)), row.names = FALSE, digits = 4)
cat("\nLimits: boundary-to-point pairs only (ladder-to-ladder pairs need the law at each point);\n",
  "the floor is a two-point lower bound on max RMSE, not an attainable rate; TV is Monte Carlo.\n", sep = "")
