#!/usr/bin/env Rscript
# End-to-end check of the SHIPPED reps > 1 analytic se -- ward_residual() as a user
# calls it -- against the exact truth of a phi_decomposition cell, with the oracle-
# nuisance se as the benchmark.
# Usage: Rscript inst/sim/se_shipped_check.R <binY 0/1> <s> <tau> [nrep=250] [reps=10]
# Pass (FINDINGS 2026-08-22, Open item 1): median(se)/empSD ~ 0.9-1.0, CV of the se
# across datasets ~0.2-0.3 (== oracle's), Wald coverage ~0.94. Output: one table,
# saved to $GAUGE_SIM_OUT/se_shipped_<cell>.rds.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
src <- readLines("inst/sim/phi_decomposition.R")
i1 <- grep("^gen <- function", src); i2 <- grep("^## ---- cells", src)
# eval(parse()) sources function DEFINITIONS from a committed repo script, not
# external input -- the standard way to reuse a CLI script's functions.
eval(parse(text = src[i1:(i2 - 1)]))

args <- commandArgs(trailingOnly = TRUE)
binY <- as.integer(args[1]) == 1L; s <- as.numeric(args[2]); tau <- as.numeric(args[3])
nrep <- if (length(args) >= 4) as.integer(args[4]) else 250L
reps <- if (length(args) >= 5) as.integer(args[5]) else 10L
n <- 800L; tr <- truth(s, tau, binY); t0 <- proc.time()[["elapsed"]]
rows <- t(vapply(seq_len(nrep), function(r) {
  seed <- 1000000L + r; set.seed(seed); d <- gen(n, s, tau, binY)
  g <- ward_residual(d, covars = "C", K = 5L, reps = reps, seed = seed * 100L, fieller = FALSE)
  so <- gauge_stats(phi_oracle(d, s, tau, binY))
  c(W = g@W, se = g@W_se, lo = g@W_ci[1], hi = g@W_ci[2],
    W_or = unname(so["W"]), se_or = unname(so["seW"]))
}, numeric(6)))
x <- as.data.frame(rows); emp <- sd(x$W); z <- qnorm(0.975)
score <- function(se, W) c(med_ratio = median(se) / emp, mean_ratio = mean(se) / emp,
  cv = sd(se) / mean(se), cov = mean(abs(W - tr$W) <= z * se),
  q05 = unname(quantile(se, .05)) / emp, q95 = unname(quantile(se, .95)) / emp)
res <- rbind(shipped = score(x$se, x$W), oracle = score(x$se_or, x$W_or))
res <- cbind(res, cov_ci = c(mean(x$lo <= tr$W & tr$W <= x$hi), NA))
cat(sprintf("shipped reps=%d analytic se: n=%d binY=%s s=%g tau=%g, nrep=%d, %.0fs | W_true=%.4f empSD=%.4f sd(W_or)=%.4f\n",
            reps, n, binY, s, tau, nrep, proc.time()[["elapsed"]] - t0, tr$W, emp, sd(x$W_or)))
print(round(res, 3))
out <- Sys.getenv("GAUGE_SIM_OUT", "gauge_sim_out"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(rows = x, res = res), file.path(out, sprintf("se_shipped_b%d_s%g_t%g_reps%d.rds", binY, s, tau, reps)))
