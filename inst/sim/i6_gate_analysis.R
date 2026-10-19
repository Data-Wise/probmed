# I6 gate (inst/sim/results/calibration/i6_authoritative_manifest_2026-10-19.csv).
# Gaussian cells (s1,s6,s7,s8): d (fresh K=200k) vs in-run exact (e): paired 90% bootstrap CI of the
# coverage difference (10000 resamples over outer reps) inside +/-0.01 and mean width ratio in [0.98, 1.02].
# Binary cells (s3,s4): report-only (b coverage vs closed-form truth, width, failed resamples).
# Usage: Rscript inst/sim/i6_gate_analysis.R <parts_dir> [expected_rows]
args <- commandArgs(TRUE); dir <- if (length(args)) args[[1]] else "inst/sim/results/calibration/i6_parts"
expected <- if (length(args) > 1) as.integer(args[[2]]) else NA_integer_
out_csv <- file.path("inst/sim/results/calibration", "i6_gate_results_2026-10-19.csv")
truth <- c(s1 = pnorm(0.4 * 0.5 / sqrt(2 * (0.16 + 1))),
           s3 = 0.5 + 0.5 * (pnorm((-2.06875 + 0.3 + 0.2) / sqrt(1.16)) - pnorm((-2.06875 + 0.3) / sqrt(1.16))),
           s4 = 0.5 + 0.5 * (pnorm((-0.25 + 0.3 + 0.2) / sqrt(1.16)) - pnorm((-0.25 + 0.3) / sqrt(1.16))),
           s6 = 0.5, s7 = 0.5, s8 = pnorm(1.5 * 2.5 / sqrt(2 * (1.5^2 + 2))))
x <- do.call(rbind, lapply(list.files(dir, pattern = "^i6_.*\\.csv$", full.names = TRUE), read.csv))
if (!is.na(expected)) stopifnot(nrow(x) == expected)
stopifnot(!anyDuplicated(x[c("scenario", "N", "rep")]), all(x$boot_attempted == 2000))
set.seed(20261019)
res <- lapply(split(x, list(x$scenario, x$N), drop = TRUE), function(cell) {
  sn <- cell$scenario[1]; N <- cell$N[1]; tp <- truth[[sn]]; n <- nrow(cell)
  cv <- function(k) as.numeric(cell[[paste0("ci_lower_", k)]] <= tp & tp <= cell[[paste0("ci_upper_", k)]])
  wd <- function(k) mean(cell[[paste0("ci_upper_", k)]] - cell[[paste0("ci_lower_", k)]])
  if (sn %in% c("s3", "s4"))
    return(data.frame(scenario = sn, N = N, role = "binary_report_only", n = n, cov_cand = mean(cv("b")),
      cov_ref = NA, diff = NA, ci_lo = NA, ci_hi = NA, width = wd("b"), width_ratio = NA,
      na_ci = sum(is.na(cell$ci_lower_b)), boot_failed = sum(cell$boot_failed_b), d_check_sd = NA,
      cov_pass = NA, width_pass = NA, pass = NA))
  dif <- cv("d") - cv("e"); bs <- replicate(10000, mean(dif[sample.int(n, n, TRUE)]))
  ci <- quantile(bs, c(.05, .95), names = FALSE); wr <- wd("d") / wd("e")
  cp <- ci[1] > -0.01 & ci[2] < 0.01; wp <- wr >= 0.98 & wr <= 1.02
  data.frame(scenario = sn, N = N, role = "gaussian_gate", n = n, cov_cand = mean(cv("d")), cov_ref = mean(cv("e")),
    diff = mean(dif), ci_lo = ci[1], ci_hi = ci[2], width = wd("d"), width_ratio = wr,
    na_ci = sum(is.na(cell$ci_lower_d)), boot_failed = sum(cell$boot_failed_d),
    d_check_sd = sd(cell$d_check_mean_z), cov_pass = cp, width_pass = wp, pass = cp & wp)
})
out <- do.call(rbind, res); rownames(out) <- NULL
write.csv(out, out_csv, row.names = FALSE)
print(out[order(out$role, out$scenario, out$N), c("scenario","N","role","n","cov_cand","cov_ref","diff","ci_lo","ci_hi","width_ratio","na_ci","boot_failed","pass")], digits = 4, row.names = FALSE)
g <- out[out$role == "gaussian_gate", ]
cat(sprintf("\nGaussian cells passing both gates: %d of %d. Live-check sd of rep-level mean z (expected 0.224): %.3f\n",
            sum(g$pass), nrow(g), mean(g$d_check_sd, na.rm = TRUE)))
cat("I6 verdict:", if (nrow(g) == 20 && all(g$pass)) "PASS (all 20 Gaussian cells)" else "FAIL or incomplete: report per-cell misses, do not alter gate", "\n")
