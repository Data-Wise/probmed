#!/usr/bin/env Rscript
# Follow-up to the frozen study's binary n=500 bias diagnostic failure.
# This script decomposes, rather than changes, the estimator or its criteria.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, dir.exists(args[1]))
out <- args[1]
devtools::load_all(quiet = TRUE)
rule <- probmed:::.wasserstein_rule(64L, normal = TRUE)
ratio <- function(p) {
  x <- sqrt(abs(p[1] - p[2]))
  y <- sqrt(abs(p[2] - p[3]))
  x / (x + y)
}
population <- function(a, ap) {
  b <- 0.8 + 0.2 * a
  mu <- -0.5 + 0.6 * a + b * 0.7 * ap
  s <- sqrt((0.5 * b + 0.3)^2 + b^2)
  stats::integrate(function(z) stats::plogis(mu + s * z) * stats::dnorm(z),
    -Inf, Inf,
    rel.tol = 1e-11
  )$value
}
target <- c(population(1, 1), population(1, 0), population(0, 0))
rows <- lapply(seq_len(2000L), function(r) {
  set.seed(20261020L + 20000L + r)
  n <- 500L
  C <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.3 * C))
  M <- 0.7 * A + 0.5 * C + stats::rnorm(n)
  eta <- -0.5 + 0.6 * A + 0.8 * M + 0.2 * A * M + 0.3 * C
  Y <- stats::rbinom(n, 1, stats::plogis(eta))
  d <- data.frame(A, M, Y, C)
  fm <- stats::lm(M ~ A + C, d)
  fy <- stats::glm(Y ~ A * M + C, d, family = stats::binomial())
  sd_m <- sqrt(mean(stats::residuals(fm)^2))
  # Independent reconstruction from fitted coefficients and conditional means.
  corner <- function(a, ap) {
    mm <- stats::predict(fm, transform(d, A = ap))
    base <- stats::predict(fy, transform(d, A = a, M = 0))
    slope <- stats::predict(fy, transform(d, A = a, M = 1)) - base
    sum(vapply(seq_along(rule$x), function(k) {
      rule$w[k] * mean(
        stats::plogis(base + slope * (mm + sd_m * rule$x[k]))
      )
    }, numeric(1)))
  }
  c(p11 = corner(1, 1), p10 = corner(1, 0), p00 = corner(0, 0))
})
ps <- do.call(rbind, rows)
utils::write.csv(ps, file.path(out, "binary-500-corner-probabilities.csv"), row.names = FALSE)
old <- utils::read.csv(file.path(out, "cell-2.csv"))
reconstructed <- apply(ps, 1L, ratio)
stopifnot(max(abs(reconstructed - old$estimate)) < 1e-8)
components <- c(
  total_bias = mean(reconstructed) - ratio(target),
  fitted_probability_bias_contribution = ratio(colMeans(ps)) - ratio(target),
  nonlinear_transformation_contribution = mean(reconstructed) - ratio(colMeans(ps)),
  reversed_indirect = sum(ps[, 1] < ps[, 2]),
  reversed_direct = sum(ps[, 2] < ps[, 3])
)
print(target)
print(components)
utils::write.csv(data.frame(component = names(components), value = unname(components)),
  file.path(out, "binary-bias-decomposition.csv"),
  row.names = FALSE
)
