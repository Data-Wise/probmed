#!/usr/bin/env Rscript
# G1d input: how often would a boundary flag "ON when |d_hat_j| / se_j < c" fail to fire (a MISS) at each
# gate-grid point? Flag OFF at a point whose near-zero contrast has noncentrality delta = d_j / se_j has probability
#   P(OFF) = P(|N(delta, 1)| >= c) = pnorm(delta - c) + pnorm(-delta - c)      [Gaussian approximation]
# At the boundary (delta = 0) this is 2 (1 - pnorm(c)). Estimator-free: uses only the DGP truths and the
# archived standard errors sqrt(v_jj) of the plug-in contrasts, not any candidate point estimator.
# Usage: Rscript scripts/tabulate-flag-misses.R <replicates.csv>   (the archived PR #55 file)
# Assumptions, stated: (1) se along each path is taken from the boundary cell (cell 5 for the direct path, v22;
# cell 6 for the indirect path, v11) at n = 500 and scaled by 1/sqrt(2) at n = 1000 (no archived boundary cell at
# n = 1000); (2) the Gaussian approximation is checked against the archive at the boundary cells 5 and 6 (control).
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, file.exists(args[1]))
d <- utils::read.csv(args[1])
d <- d[d$method == "k4_fieller", ]
cs <- c(1, 1.645, 1.96, 2.576, 3)

corner_mean <- function(direct, alpha, inter, a, ap) {
  b <- .8 + inter * a
  mu <- -.5 + direct * a + b * alpha * ap
  sdv <- sqrt((.5 * b + .3)^2 + b^2)
  stats::integrate(function(z) stats::plogis(mu + sdv * z) * stats::dnorm(z), -Inf, Inf, rel.tol = 1e-11)$value
}
d2_direct <- function(dir) corner_mean(dir, .7, 0, 1, 0) - corner_mean(dir, .7, 0, 0, 0)
d1_alpha <- function(al) corner_mean(.6, al, .2, 1, 1) - corner_mean(.6, al, .2, 1, 0)

se_of <- function(cell, v) sqrt(stats::median(d[d$cell == cell, v], na.rm = TRUE))
se2 <- se_of(5, "v22"); se1 <- se_of(6, "v11")
# constancy check of the se along the path: cell 7 (direct .05) vs cell 5 (direct 0); cell 2 (alpha .7) vs cell 6 (alpha 0)
se2_c7 <- se_of(7, "v22"); se1_c2 <- se_of(2, "v11")
cat(sprintf("Archived SE at n = 500: direct path se(d2) %.4f (boundary cell 5), %.4f at cell 7 (direct .05); indirect path se(d1) %.4f (boundary cell 6), %.4f at cell 2 (alpha .7)\n",
  se2, se2_c7, se1, se1_c2))

# control: Gaussian P(OFF) at delta = 0 versus the archived boundary cells
emp <- function(cell, dj, vj) { x <- abs(d[d$cell == cell, dj]) / sqrt(d[d$cell == cell, vj]); sapply(cs, function(c) mean(x >= c, na.rm = TRUE)) }
ctl <- rbind(c = cs, gaussian = 2 * (1 - stats::pnorm(cs)), cell5_d2 = emp(5, "d2", "v22"), cell6_d1 = emp(6, "d1", "v11"))
cat("\nControl at the boundary (true contrast 0): P(OFF) = P(|d_hat|/se >= c); Gaussian vs archived (2000 replicates)\n")
print(round(ctl, 4))

grids <- list(
  `500`  = list(direct = c(0, .00875, .0175, .02635, .0352, .05325, .0713, .1095, .1477),
                indirect = c(0, .0040, .0080, .0120, .0160, .0242, .0324, .05005, .0677)),
  `1000` = list(direct = c(0, .00625, .0125, .01875, .025, .03775, .0505, .0777, .1049),
                indirect = c(0, .00285, .0057, .00855, .0114, .01725, .0231, .03525, .0474)))
poff <- function(delta, c) stats::pnorm(delta - c) + stats::pnorm(-delta - c)
for (n in c(500, 1000)) {
  sc <- if (n == 1000) sqrt(.5) else 1
  for (path in c("direct", "indirect")) {
    g <- grids[[as.character(n)]][[path]]
    dj <- if (path == "direct") vapply(g, d2_direct, 0) else vapply(g, d1_alpha, 0)
    delta <- abs(dj) / ((if (path == "direct") se2 else se1) * sc)
    tab <- data.frame(param = g, true_contrast = dj, delta = delta)
    for (c in cs) tab[[sprintf("OFF@%.3g", c)]] <- poff(delta, c)
    cat(sprintf("\nn = %d, %s path: P(flag OFF) = miss probability by threshold c (SE units)\n", n, path))
    print(tab, row.names = FALSE, digits = 3)
  }
}
cat("\nReading: a MISS is the flag staying OFF while the near-zero contrast is small relative to its SE; at delta = 0 it is the\n",
  "size 2 (1 - pnorm(c)). Whether a point should be flagged is the author's definition (G1d); this table only gives probabilities.\n", sep = "")
