#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, dir.exists(args[1]))
out <- args[1]
devtools::load_all(quiet = TRUE)
source("scripts/wasserstein-boundary-candidates.R")
started <- Sys.time()
cells <- data.frame(
  cell = 1:7, binary = c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
  n = c(500L, 500L, 1000L, 1000L, 500L, 500L, 500L),
  scenario = c(rep("regular", 4), "zero_direct", "zero_indirect", "weak_direct"),
  alpha = c(rep(0.7, 5), 0, 0.7),
  direct = c(rep(0.6, 4), 0, 0.6, 0.05),
  interaction = c(rep(0.2, 4), 0, 0.2, 0)
)
truth <- function(s) {
  corner <- function(a, ap) {
    b <- 0.8 + s$interaction * a
    mu <- -0.5 + s$direct * a + b * s$alpha * ap
    sd <- sqrt((0.5 * b + 0.3)^2 + b^2 + !s$binary)
    if (s$binary) {
      stats::integrate(function(z) stats::plogis(mu + sd * z) * stats::dnorm(z),
        -Inf, Inf,
        rel.tol = 1e-11
      )$value
    } else {
      c(mu, sd)
    }
  }
  cs <- list(corner(1, 1), corner(1, 0), corner(0, 0))
  if (s$binary) {
    return(wb_ratio(c(cs[[1]] - cs[[2]], cs[[2]] - cs[[3]])))
  }
  w <- c(sqrt(sum((cs[[1]] - cs[[2]])^2)), sqrt(sum((cs[[2]] - cs[[3]])^2)))
  w[1] / sum(w)
}
run <- function(cell) {
  s <- cells[cell, ]
  target <- truth(s)
  methods <- if (s$binary) c("raw", "damped_k4", "damped_k9") else "raw"
  rows <- lapply(seq_len(2000L), function(r) {
    seed <- 20361020L + cell * 10000L + r
    set.seed(seed)
    C <- stats::rnorm(s$n)
    A <- stats::rbinom(s$n, 1, stats::plogis(0.3 * C))
    M <- s$alpha * A + 0.5 * C + stats::rnorm(s$n)
    eta <- -0.5 + s$direct * A + 0.8 * M + s$interaction * A * M + 0.3 * C
    Y <- if (s$binary) stats::rbinom(s$n, 1, stats::plogis(eta)) else eta + stats::rnorm(s$n)
    d <- data.frame(A, M, Y, C)
    result <- tryCatch(
      {
        if (s$binary) {
          fit <- wb_binary(d)
          if (r <= 10L) {
            check <- probmed:::.wasserstein_observed(d, binY = TRUE)
            stopifnot(abs(check$pmedW - fit$raw) < 1e-9, abs(check$se - fit$se) < 1e-8)
          }
          estimates <- c(fit$raw, wb_correct(fit$d, fit$v, 4), wb_correct(fit$d, fit$v, 9))
        } else {
          fit <- probmed:::.wasserstein_observed(d)
          fit$min_z <- NA_real_
          estimates <- fit$pmedW
        }
        data.frame(
          method = methods, estimate = estimates, se = fit$se,
          min_z = fit$min_z, error = ""
        )
      },
      error = function(e) {
        data.frame(
          method = methods, estimate = NA_real_, se = NA_real_,
          min_z = NA_real_, error = conditionMessage(e)
        )
      }
    )
    cbind(cell,
      scenario = s$scenario, binary = s$binary, n = s$n,
      replicate = r, seed, truth = target, result
    )
  })
  result <- do.call(rbind, rows)
  utils::write.csv(result, file.path(out, paste0("cell-", cell, ".csv")), row.names = FALSE)
  result
}
results <- do.call(rbind, parallel::mclapply(seq_len(7L), run, mc.cores = 4L))
utils::write.csv(results, file.path(out, "replicates.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(results, interaction(results$cell, results$method, drop = TRUE)), function(d) {
  good <- is.finite(d$estimate) & is.finite(d$se) & d$error == ""
  x <- d[good, ]
  sd <- stats::sd(x$estimate)
  bias <- mean(x$estimate - x$truth)
  coverage <- mean(abs(x$estimate - x$truth) <= stats::qnorm(0.975) * x$se)
  ratio <- mean(x$se) / sd
  data.frame(
    cell = d$cell[1], scenario = d$scenario[1], n = d$n[1],
    method = d$method[1], reps = nrow(d), failed = sum(!good), truth = d$truth[1],
    bias, bias_mcse = sd / sqrt(nrow(x)), empirical_sd = sd,
    standardized_bias = abs(bias) / sd, mean_se = mean(x$se), se_ratio = ratio,
    coverage, coverage_mcse = sqrt(coverage * (1 - coverage) / nrow(x)),
    rmse = sqrt(mean((x$estimate - x$truth)^2)),
    outside = sum(x$estimate < 0 | x$estimate > 1),
    near_zero = if (all(is.na(x$min_z))) NA_real_ else mean(x$min_z < 2),
    pass = all(good) && abs(bias) / sd < 0.10 && ratio >= 0.90 && ratio <= 1.10 &&
      coverage >= 0.92 && coverage <= 0.98
  )
}))
utils::write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
writeLines(c(
  capture.output(sessionInfo()),
  paste("Elapsed seconds:", as.numeric(difftime(Sys.time(), started, units = "secs"))),
  paste("Source HEAD:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
  paste("Candidate SHA256:", system2("shasum",
    c("-a", "256", "scripts/wasserstein-boundary-candidates.R"), stdout = TRUE
  ))
), file.path(out, "provenance.txt"))
print(summary)
if (!all(summary$pass)) message("Some frozen diagnostics fail; merge remains blocked.")
