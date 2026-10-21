#!/usr/bin/env Rscript
# Targeted post-pilot diagnostic, not a validation sample.
devtools::load_all(quiet = TRUE)
set.seed(20261020L + 20000L + 166L)
n <- 500L
C <- stats::rnorm(n)
A <- stats::rbinom(n, 1, stats::plogis(0.3 * C))
M <- 0.7 * A + 0.5 * C + stats::rnorm(n)
eta <- -0.5 + 0.6 * A + 0.8 * M + 0.2 * A * M + 0.3 * C
Y <- stats::rbinom(n, 1, stats::plogis(eta))
d <- data.frame(A, M, Y, C)
fit <- function(d) {
  probmed:::.wasserstein_observed(d,
    binY = TRUE, influence = FALSE
  )$pmedW
}
full <- fit(d)
deleted <- vapply(seq_len(n), function(i) fit(d[-i, ]), numeric(1))
print(c(
  replicate = 166, raw = full, delete_one = n * full - (n - 1) * mean(deleted),
  minimum_deleted = min(deleted), maximum_deleted = max(deleted)
))
