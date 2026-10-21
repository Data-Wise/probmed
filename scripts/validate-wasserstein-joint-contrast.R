#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, dir.exists(args[1]))
out <- args[1]
devtools::load_all(quiet = TRUE)
source("scripts/wasserstein-joint-contrast.R")
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
    delta <- c(cs[[1]] - cs[[2]], cs[[2]] - cs[[3]])
    return(list(p = wb_ratio(delta), d = delta))
  }
  w <- c(sqrt(sum((cs[[1]] - cs[[2]])^2)), sqrt(sum((cs[[2]] - cs[[3]])^2)))
  list(p = w[1] / sum(w), d = c(NA_real_, NA_real_))
}
run <- function(cell) {
  s <- cells[cell, ]
  target <- truth(s)
  methods <- if (s$binary) c("raw_wald", "raw_projected", "k4_projected") else "raw_wald"
  rows <- lapply(seq_len(2000L), function(r) {
    seed <- 20461020L + cell * 10000L + r
    set.seed(seed)
    C <- stats::rnorm(s$n)
    A <- stats::rbinom(s$n, 1, stats::plogis(0.3 * C))
    M <- s$alpha * A + 0.5 * C + stats::rnorm(s$n)
    eta <- -0.5 + s$direct * A + 0.8 * M + s$interaction * A * M + 0.3 * C
    Y <- if (s$binary) stats::rbinom(s$n, 1, stats::plogis(eta)) else eta + stats::rnorm(s$n)
    result <- tryCatch(
      {
        data <- data.frame(A, M, Y, C)
        if (s$binary) {
          fit <- wb_binary(data)
          if (r <= 10L) {
            check <- probmed:::.wasserstein_observed(data, binY = TRUE)
            stopifnot(abs(check$pmedW - fit$raw) < 1e-9, abs(check$se - fit$se) < 1e-8)
          }
          ci <- wj_interval(fit$d, fit$v)
          estimates <- c(fit$raw, fit$raw, wb_correct(fit$d, fit$v, 4))
          lower <- c(fit$raw - stats::qnorm(0.975) * fit$se, rep(ci[1], 2))
          upper <- c(fit$raw + stats::qnorm(0.975) * fit$se, rep(ci[2], 2))
          difference <- fit$d - target$d
          joint <- drop(crossprod(difference, solve(fit$v, difference))) <= stats::qchisq(0.95, 2)
        } else {
          fit <- probmed:::.wasserstein_observed(data)
          estimates <- fit$pmedW
          lower <- estimates - stats::qnorm(0.975) * fit$se
          upper <- estimates + stats::qnorm(0.975) * fit$se
          fit$d <- rep(NA_real_, 2)
          fit$v <- matrix(NA_real_, 2, 2)
          joint <- NA
        }
        data.frame(
          method = methods, estimate = estimates, se = fit$se, lower, upper,
          d1 = fit$d[1], d2 = fit$d[2], v11 = fit$v[1, 1], v12 = fit$v[1, 2],
          v22 = fit$v[2, 2], joint_covered = joint, error = ""
        )
      },
      error = function(e) {
        data.frame(
          method = methods, estimate = NA_real_, se = NA_real_, lower = NA_real_,
          upper = NA_real_, d1 = NA_real_, d2 = NA_real_, v11 = NA_real_, v12 = NA_real_,
          v22 = NA_real_, joint_covered = NA, error = conditionMessage(e)
        )
      }
    )
    cbind(cell,
      scenario = s$scenario, n = s$n, replicate = r, seed,
      truth = target$p, true_d1 = target$d[1], true_d2 = target$d[2], result
    )
  })
  result <- do.call(rbind, rows)
  utils::write.csv(result, file.path(out, paste0("cell-", cell, ".csv")), row.names = FALSE)
  message("Completed cell ", cell)
  result
}
results <- do.call(rbind, parallel::mclapply(seq_len(7L), run, mc.cores = 4L))
utils::write.csv(results, file.path(out, "replicates.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(results, interaction(results$cell, results$method, drop = TRUE)), function(d) {
  good <- is.finite(d$estimate) & is.finite(d$se) & is.finite(d$lower) & is.finite(d$upper) & d$error == ""
  x <- d[good, ]
  sd <- stats::sd(x$estimate)
  bias <- mean(x$estimate - x$truth)
  coverage <- mean(x$lower <= x$truth & x$upper >= x$truth)
  ratio <- mean(x$se) / sd
  data.frame(
    cell = d$cell[1], scenario = d$scenario[1], method = d$method[1], reps = nrow(d),
    failed = sum(!good), bias, bias_mcse = sd / sqrt(nrow(x)), empirical_sd = sd,
    standardized_bias = abs(bias) / sd, se_ratio = ratio, coverage,
    coverage_mcse = sqrt(coverage * (1 - coverage) / nrow(x)),
    mean_width = mean(x$upper - x$lower), full_range = mean(x$lower == 0 & x$upper == 1),
    joint_coverage = if (all(is.na(x$joint_covered))) NA_real_ else mean(x$joint_covered),
    point_in_interval = mean(x$lower <= x$estimate & x$upper >= x$estimate),
    pass = all(good) && abs(bias) / sd < 0.10 && ratio >= 0.90 && ratio <= 1.10 &&
      coverage >= 0.92 && coverage <= 0.98
  )
}))
utils::write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
writeLines(c(
  capture.output(sessionInfo()),
  paste("Elapsed seconds:", as.numeric(difftime(Sys.time(), started, units = "secs"))),
  paste("Source HEAD:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
  system2("shasum", c(
    "-a", "256", "scripts/wasserstein-joint-contrast.R",
    "scripts/wasserstein-boundary-candidates.R", "scripts/validate-wasserstein-joint-contrast.R"
  ), stdout = TRUE)
), file.path(out, "provenance.txt"))
print(summary)
if (!all(summary$pass)) message("Frozen diagnostics fail; both PRs remain blocked.")
