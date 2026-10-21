#!/usr/bin/env Rscript
# Post-weight-fix re-measurement of ward_residual()'s bootstrap arm and the two gates
# (weak_id = A2, oe_regular = A1) against exact truth. One cell per call.
# Run from the package root. Usage: Rscript inst/sim/boot_gates_check.R <binY 0/1> <s> <tau> [nrep=200] [B=200]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
src <- readLines("inst/sim/phi_decomposition.R")
i1 <- grep("^gen <- function", src); i2 <- grep("^## ---- cells", src)
# eval(parse()) sources function DEFINITIONS (gen/truth/oracle) from a committed repo script.
eval(parse(text = src[i1:(i2 - 1)]))
args <- commandArgs(trailingOnly = TRUE)
binY <- as.integer(args[1]) == 1L; s <- as.numeric(args[2]); tau <- as.numeric(args[3])
nrep <- if (length(args) >= 4) as.integer(args[4]) else 200L
B    <- if (length(args) >= 5) as.integer(args[5]) else 200L
n <- 800L; tr <- truth(s, tau, binY); t0 <- proc.time()[["elapsed"]]
rows <- t(vapply(seq_len(nrep), function(r) {
  seed <- 1000000L + r; set.seed(seed); d <- gen(n, s, tau, binY)
  g <- suppressWarnings(ward_residual(d, covars = "C", K = 5L, se_method = "bootstrap", B = B,
                                      seed = seed * 100L, fieller = TRUE))
  so <- gauge_stats(phi_oracle(d, s, tau, binY))
  if (r %% 25 == 0) message(sprintf("  %d/%d  %.0fs", r, nrep, proc.time()[["elapsed"]] - t0))
  c(W = g@W, se_boot = g@W_se, lo_pct = g@W_ci[1], hi_pct = g@W_ci[2],
    lo_wald = g@W_ci_wald[1], hi_wald = g@W_ci_wald[2],
    weak_id = as.numeric(g@weak_id), wid_ratio = g@weak_id_ratio,
    oe_snr = g@oe_snr, oe_regular = as.numeric(g@oe_regular),
    fieller_bounded = as.numeric(identical(g@fieller_type, "bounded")),
    W_or = unname(so["W"]), se_or = unname(so["seW"]))
}, numeric(13)))
x <- as.data.frame(rows); emp <- sd(x$W); Wt <- tr$W
inside <- function(lo, hi) lo <= Wt & Wt <= hi
cov_pct <- inside(x$lo_pct, x$hi_pct); cov_wald <- inside(x$lo_wald, x$hi_wald)
se_wald <- (x$hi_wald - x$lo_wald) / (2 * qnorm(0.975))
wr <- (x$hi_pct - x$lo_pct) / (x$hi_wald - x$lo_wald)
flag <- x$weak_id == 1; reg <- x$oe_regular == 1
cat(sprintf("cell: n=%d binY=%s s=%g tau=%g | nrep=%d B=%d %.0fs | W_true=%.4f empSD=%.4f sd(W_or)=%.4f\n",
            n, binY, s, tau, nrep, B, proc.time()[["elapsed"]] - t0, Wt, emp, sd(x$W_or)))
out <- c(cov_percentile = mean(cov_pct), cov_wald = mean(cov_wald),
         oracle_cov = mean(abs(x$W_or - Wt) <= qnorm(0.975) * x$se_or),
         boot_sd_over_empSD_med = median(x$se_boot) / emp, boot_sd_cv = sd(x$se_boot) / mean(x$se_boot),
         wald_se_over_empSD_med = median(se_wald) / emp,
         pct_width_over_calib_med = median(x$hi_pct - x$lo_pct) / (2 * qnorm(0.975) * emp),
         width_ratio_pct_wald_med = median(wr), width_ratio_q90 = unname(quantile(wr, .9)),
         weak_id_rate = mean(flag), oe_regular_rate = mean(reg), oe_snr_med = median(x$oe_snr),
         fieller_bounded_rate = mean(x$fieller_bounded == 1),
         cov_wald_flagged = if (any(flag)) mean(cov_wald[flag]) else NA,
         cov_wald_unflagged = if (any(!flag)) mean(cov_wald[!flag]) else NA,
         cov_pct_flagged = if (any(flag)) mean(cov_pct[flag]) else NA,
         cov_wald_oe_irregular = if (any(!reg)) mean(cov_wald[!reg]) else NA,
         cov_wald_oe_regular = if (any(reg)) mean(cov_wald[reg]) else NA,
         n_flagged = sum(flag), n_irregular = sum(!reg))
print(round(out, 3))
out_dir <- Sys.getenv("GAUGE_SIM_OUT", "gauge_sim_out"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(rows = x, summary = out, cell = list(n = n, binY = binY, s = s, tau = tau, nrep = nrep, B = B)),
        file.path(out_dir, sprintf("boot_gates_b%d_s%g_t%g.rds", binY, s, tau)))
