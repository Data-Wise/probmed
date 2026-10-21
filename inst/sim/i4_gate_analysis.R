# I4 prespecified gate (tasks/plan-mc-bootstrap.md, I4): paired coverage equivalence
# (90% paired bootstrap CI over outer reps, 10000 resamples, inside +/-0.01) and
# mean-width ratio in [0.98, 1.02], per cell, vs exact (Gaussian) or d at K=1e6 (binary).
dir <- "inst/sim/results/calibration"
# Closed-form true P_med (independent cross-world arms; see SCENARIOS in the pmed sim config).
# truth.csv holds 2e6-draw MC values (SE ~3.5e-4), too coarse for the null cells (exactly 0.5).
truth <- c(s1 = pnorm(0.4 * 0.5 / sqrt(2 * (0.16 + 1))),
           s3 = 0.5 + 0.5 * (pnorm((-2.06875 + 0.3 + 0.2) / sqrt(1.16)) - pnorm((-2.06875 + 0.3) / sqrt(1.16))),
           s4 = 0.5 + 0.5 * (pnorm((-0.25 + 0.3 + 0.2) / sqrt(1.16)) - pnorm((-0.25 + 0.3) / sqrt(1.16))),
           s6 = 0.5, s7 = 0.5,
           s8 = pnorm(1.5 * 2.5 / sqrt(2 * (1.5^2 + 2))))
x <- do.call(rbind, lapply(list.files(file.path(dir, "i4_parts"), full.names = TRUE), read.csv))
stopifnot(nrow(x) == 180 * 200, !anyDuplicated(x[c("scenario", "N", "rep")]))
set.seed(20261019)
res <- list()
for (cell in split(x, list(x$scenario, x$N), drop = TRUE)) {
  sn <- cell$scenario[1]; N <- cell$N[1]; tp <- truth[[sn]]
  bin <- sn %in% c("s3", "s4"); ref <- if (bin) "d" else "e"
  lo <- function(k) cell[[paste0("ci_lower_", k)]]; hi <- function(k) cell[[paste0("ci_upper_", k)]]
  cov <- function(k) as.numeric(lo(k) <= tp & tp <= hi(k))
  wid <- function(k) hi(k) - lo(k)
  for (k in setdiff(c("a", "b", "c", "d"), ref)) {
    dif <- cov(k) - cov(ref); n <- length(dif)
    bs <- replicate(10000, mean(dif[sample.int(n, n, TRUE)]))
    ci <- quantile(bs, c(.05, .95), names = FALSE)
    wr <- mean(wid(k), na.rm = TRUE) / mean(wid(ref), na.rm = TRUE)
    res[[length(res) + 1]] <- data.frame(scenario = sn, N = N, cand = k, ref = ref,
      n = n, na_ci = sum(is.na(lo(k))), cov_cand = mean(cov(k)), cov_ref = mean(cov(ref)),
      diff = mean(dif), ci_lo = ci[1], ci_hi = ci[2], width_ratio = wr,
      cov_pass = ci[1] > -0.01 & ci[2] < 0.01, width_pass = wr >= 0.98 & wr <= 1.02,
      boot_failed = sum(cell[[paste0("boot_failed_", if (k == "d") "d" else k)]] %||% 0, na.rm = TRUE))
  }
}
`%||%` <- function(a, b) if (is.null(a)) b else a
out <- do.call(rbind, res); out$pass <- out$cov_pass & out$width_pass
write.csv(out, file.path(dir, "i4_gate_results_2026-10-19.csv"), row.names = FALSE)
print(out[order(out$cand, out$scenario, out$N), c(1:3, 7:14)], digits = 4, row.names = FALSE)
cat("\nPass in ALL cells, by candidate:\n"); print(tapply(out$pass, out$cand, all))
cat("Cells passed per candidate:\n"); print(tapply(out$pass, out$cand, sum))
