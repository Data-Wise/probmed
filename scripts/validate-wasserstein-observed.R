#!/usr/bin/env Rscript
# Prespecified design: docs/reports/PLAN-wasserstein-observed-inference-2026-10-20.md
# Run from the package root. Results are external to the clean CI checkout.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else stop("Supply an existing output directory")
stopifnot(dir.exists(out))
devtools::load_all(quiet = TRUE)
started <- Sys.time()
cells <- expand.grid(binary = c(FALSE, TRUE), n = c(500L, 1000L))
truth <- function(binary) {
  corner <- function(a, ap) {
    b <- 0.8 + 0.2 * a
    mu <- -0.5 + 0.6 * a + b * 0.7 * ap
    s <- sqrt((0.5 * b + 0.3)^2 + b^2 + if (binary) 0 else 1)
    if (binary) {
      stats::integrate(function(z) {
        stats::plogis(mu + s * z) *
          stats::dnorm(z)
      }, -Inf, Inf, rel.tol = 1e-11)$value
    } else {
      c(mu = mu, sd = s)
    }
  }
  cs <- list(corner(1, 1), corner(1, 0), corner(0, 0))
  dist <- function(i, j) {
    if (binary) {
      sqrt(abs(cs[[i]] - cs[[j]]))
    } else {
      sqrt(sum((cs[[i]] - cs[[j]])^2))
    }
  }
  dist(1, 2) / (dist(1, 2) + dist(2, 3))
}
run <- function(cell) {
  binary <- cells$binary[cell]
  n <- cells$n[cell]
  target <- truth(binary)
  rows <- lapply(seq_len(2000L), function(r) {
    seed <- 20261020L + cell * 10000L + r
    set.seed(seed)
    C <- stats::rnorm(n)
    A <- stats::rbinom(n, 1, stats::plogis(0.3 * C))
    M <- 0.7 * A + 0.5 * C + stats::rnorm(n)
    eta <- -0.5 + 0.6 * A + 0.8 * M + 0.2 * A * M + 0.3 * C
    Y <- if (binary) stats::rbinom(n, 1, stats::plogis(eta)) else eta + stats::rnorm(n)
    fit <- tryCatch(
      probmed:::.wasserstein_observed(data.frame(A, M, Y, C), binY = binary),
      error = identity
    )
    ok <- !inherits(fit, "error")
    data.frame(cell, binary, n,
      replicate = r, seed, truth = target,
      estimate = if (ok) fit$pmedW else NA_real_,
      se = if (ok) fit$se else NA_real_,
      error = if (ok) "" else conditionMessage(fit)
    )
  })
  result <- do.call(rbind, rows)
  utils::write.csv(result, file.path(out, paste0("cell-", cell, ".csv")), row.names = FALSE)
  result
}
# Separate processes only read the frozen source; do not edit while running.
results <- do.call(rbind, parallel::mclapply(seq_len(nrow(cells)), run, mc.cores = 4L))
utils::write.csv(results, file.path(out, "replicates.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(results, results$cell), function(d) {
  good <- is.finite(d$estimate) & is.finite(d$se)
  x <- d[good, ]
  sd <- stats::sd(x$estimate)
  bias <- mean(x$estimate) - x$truth[1]
  coverage <- mean(abs(x$estimate - x$truth) <= stats::qnorm(0.975) * x$se)
  data.frame(
    binary = d$binary[1], n = d$n[1], reps = nrow(d), failed = sum(!good),
    truth = d$truth[1], bias, bias_mcse = sd / sqrt(nrow(x)), empirical_sd = sd,
    mean_se = mean(x$se), se_ratio = mean(x$se) / sd, coverage,
    coverage_mcse = sqrt(coverage * (1 - coverage) / nrow(x)),
    pass = all(good) && abs(bias) / sd < 0.10 && mean(x$se) / sd >= 0.90 &&
      mean(x$se) / sd <= 1.10 && coverage >= 0.92 && coverage <= 0.98
  )
}))
utils::write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
writeLines(
  c(
    capture.output(sessionInfo()), paste(
      "Elapsed seconds:",
      as.numeric(difftime(Sys.time(), started, units = "secs"))
    ),
    paste("Source HEAD:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
  ),
  file.path(out, "provenance.txt")
)
print(summary)
if (!all(summary$pass)) stop("Prespecified diagnostics failed; investigate without changing gates")
