#!/usr/bin/env Rscript
# Exploratory only: never modifies the public estimator or validation gates.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, dir.exists(args[1]))
devtools::load_all(quiet = TRUE)
truth <- utils::read.csv(file.path(args[1], "summary.csv"))
truth <- truth$truth[truth$binary & truth$n == 500L]
stopifnot(length(truth) == 1L)
rows <- lapply(seq_len(200L), function(r) {
  set.seed(20261020L + 20000L + r)
  n <- 500L
  C <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.3 * C))
  M <- 0.7 * A + 0.5 * C + stats::rnorm(n)
  eta <- -0.5 + 0.6 * A + 0.8 * M + 0.2 * A * M + 0.3 * C
  Y <- stats::rbinom(n, 1, stats::plogis(eta))
  d <- data.frame(A, M, Y, C)
  tryCatch(
    {
      full <- probmed:::.wasserstein_observed(d, binY = TRUE)
      group <- rep(seq_len(20L), length.out = n)
      deleted <- vapply(seq_len(20L), function(g) {
        probmed:::.wasserstein_observed(d[group != g, ],
          binY = TRUE, influence = FALSE
        )$pmedW
      }, numeric(1))
      data.frame(
        replicate = r, raw = full$pmedW,
        corrected = 20 * full$pmedW - 19 * mean(deleted),
        se = full$se, error = ""
      )
    },
    error = function(e) {
      data.frame(
        replicate = r, raw = NA_real_, corrected = NA_real_,
        se = NA_real_, error = conditionMessage(e)
      )
    }
  )
})
results <- do.call(rbind, rows)
utils::write.csv(results, file.path(args[1], "bias-correction-pilot.csv"),
  row.names = FALSE
)
stopifnot(all(results$error == ""))
summarize <- function(x) {
  c(
    bias = mean(x - truth), bias_mcse = stats::sd(x) / sqrt(length(x)),
    sd = stats::sd(x), rmse = sqrt(mean((x - truth)^2)),
    outside = sum(x < 0 | x > 1),
    exploratory_coverage = mean(abs(x - truth) <= stats::qnorm(0.975) * results$se)
  )
}
summary <- rbind(raw = summarize(results$raw), corrected = summarize(results$corrected))
print(summary)
print(c(
  mean_correction = mean(results$corrected - results$raw),
  correction_sd = stats::sd(results$corrected - results$raw),
  correction_mcse = stats::sd(results$corrected - results$raw) / sqrt(nrow(results))
))
utils::write.csv(summary, file.path(args[1], "bias-correction-pilot-summary.csv"))
