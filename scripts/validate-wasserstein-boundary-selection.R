#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, dir.exists(args[1]))
out <- args[1]
devtools::load_all(quiet = TRUE)
source("scripts/wasserstein-boundary-selection.R")
# Reuse only the frozen DGP and truth definitions, never its executing driver.
previous <- parse("scripts/validate-wasserstein-joint-contrast.R")
for (name in c("cells", "truth")) {
  definition <- Filter(function(e) {
    is.call(e) && identical(e[[1]], as.name("<-")) &&
      identical(e[[2]], as.name(name))
  }, as.list(previous))
  stopifnot(length(definition) == 1L)
  eval(definition[[1]])
}
set.seed(816)
z <- matrix(stats::rnorm(32768 * 2), ncol = 2)
z <- rbind(z, -z)
started <- Sys.time()
run <- function(cell) {
  s <- cells[cell, ]
  target <- truth(s)
  methods <- if (s$binary) c("raw_wald", "k4_fieller", "select95_fieller", "select995_fieller") else "raw_wald"
  result <- do.call(rbind, lapply(seq_len(2000L), function(r) {
    seed <- 20561020L + cell * 10000L + r
    set.seed(seed)
    C <- stats::rnorm(s$n)
    A <- stats::rbinom(s$n, 1, stats::plogis(.3 * C))
    M <- s$alpha * A + .5 * C + stats::rnorm(s$n)
    eta <- -.5 + s$direct * A + .8 * M + s$interaction * A * M + .3 * C
    Y <- if (s$binary) stats::rbinom(s$n, 1, stats::plogis(eta)) else eta + stats::rnorm(s$n)
    value <- tryCatch(
      {
        data <- data.frame(A, M, Y, C)
        if (s$binary) {
          fit <- wb_binary(data)
          if (r <= 10L) {
            oracle <- probmed:::.wasserstein_observed(data, binY = TRUE)
            stopifnot(abs(oracle$pmedW - fit$raw) < 1e-9, abs(oracle$se - fit$se) < 1e-8)
          }
          ci <- ws_interval(fit$d, fit$v)
          a <- ws_fit(fit$d, fit$v, stats::qnorm(.975), z)
          b <- ws_fit(fit$d, fit$v, stats::qnorm(.9975), z)
          estimates <- c(fit$raw, wb_correct(fit$d, fit$v), a$estimate, b$estimate)
          ses <- c(fit$se, fit$se, a$se, b$se)
          lower <- c(fit$raw - stats::qnorm(.975) * fit$se, rep(ci[1], 3))
          upper <- c(fit$raw + stats::qnorm(.975) * fit$se, rep(ci[2], 3))
          selected <- c(0L, 0L, a$selected, b$selected)
          ambiguous <- c(FALSE, FALSE, a$ambiguous, b$ambiguous)
          center1 <- c(fit$d[1], fit$d[1], a$center[1], b$center[1])
          center2 <- c(fit$d[2], fit$d[2], a$center[2], b$center[2])
        } else {
          fit <- probmed:::.wasserstein_observed(data)
          estimates <- fit$pmedW
          ses <- fit$se
          lower <- estimates - stats::qnorm(.975) * ses
          upper <- estimates + stats::qnorm(.975) * ses
          fit$d <- rep(NA_real_, 2)
          fit$v <- matrix(NA_real_, 2, 2)
          selected <- 0L
          ambiguous <- FALSE
          center1 <- center2 <- NA_real_
        }
        data.frame(
          method = methods, estimate = estimates, se = ses, lower, upper,
          d1 = fit$d[1], d2 = fit$d[2], v11 = fit$v[1, 1], v12 = fit$v[1, 2],
          v22 = fit$v[2, 2], selected, ambiguous, center1, center2, error = ""
        )
      },
      error = function(e) {
        data.frame(
          method = methods, estimate = NA_real_, se = NA_real_, lower = NA_real_,
          upper = NA_real_, d1 = NA_real_, d2 = NA_real_, v11 = NA_real_, v12 = NA_real_,
          v22 = NA_real_, selected = NA_integer_, ambiguous = NA, center1 = NA_real_,
          center2 = NA_real_, error = conditionMessage(e)
        )
      }
    )
    cbind(cell,
      scenario = s$scenario, n = s$n, replicate = r, seed,
      truth = target$p, true_d1 = target$d[1], true_d2 = target$d[2], value
    )
  }))
  utils::write.csv(result, file.path(out, paste0("cell-", cell, ".csv")), row.names = FALSE)
  message("Completed cell ", cell)
  result
}
results <- do.call(rbind, parallel::mclapply(1:7, run, mc.cores = 4L))
utils::write.csv(results, file.path(out, "replicates.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(results, interaction(results$cell, results$method, drop = TRUE)), function(d) {
  good <- is.finite(d$estimate) & is.finite(d$se) & is.finite(d$lower) & is.finite(d$upper) & d$error == ""
  x <- d[good, ]
  sd <- stats::sd(x$estimate)
  bias <- mean(x$estimate - x$truth)
  coverage <- mean(x$lower <= x$truth & x$upper >= x$truth)
  ratio <- mean(x$se) / sd
  data.frame(
    cell = d$cell[1], method = d$method[1], reps = nrow(d), failed = sum(!good),
    bias, bias_mcse = sd / sqrt(nrow(x)), empirical_sd = sd, standardized_bias = abs(bias) / sd,
    se_ratio = ratio, coverage, coverage_mcse = sqrt(coverage * (1 - coverage) / nrow(x)),
    mean_width = mean(x$upper - x$lower), selected_rate = mean(x$selected != 0),
    ambiguous_rate = mean(x$ambiguous), point_in_interval = mean(x$lower <= x$estimate & x$upper >= x$estimate),
    pass = all(good) && is.finite(sd) && sd > 0 && abs(bias) / sd < .10 &&
      ratio >= .90 && ratio <= 1.10 && coverage >= .92 && coverage <= .98
  )
}))
utils::write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
writeLines(
  c(
    capture.output(sessionInfo()), paste("Elapsed seconds:", as.numeric(difftime(Sys.time(), started, units = "secs"))),
    paste("Source HEAD:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    system2("shasum", c(
      "-a", "256", "scripts/wasserstein-boundary-selection.R",
      "scripts/wasserstein-boundary-candidates.R", "scripts/wasserstein-joint-contrast.R",
      "scripts/validate-wasserstein-boundary-selection.R", "scripts/validate-wasserstein-joint-contrast.R"
    ), stdout = TRUE)
  ),
  file.path(out, "provenance.txt")
)
print(summary)
if (!all(summary$pass)) message("Frozen diagnostics fail; both PRs remain blocked.")
