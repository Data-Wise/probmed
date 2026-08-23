# Calibration analysis over the boot_gates_check.R outputs. Usage: Rscript inst/sim/boot_gates_calibrate.R <dir with boot_gates_*.rds>
S <- commandArgs(trailingOnly = TRUE)[1]
cells <- list(cont_strong = "boot_gates_b0_s1_t0.rds", bin_strong = "boot_gates_b1_s1_t0.rds",
              intermediate = "boot_gates_b0_s0.5_t0.4.rds", near_null = "boot_gates_b0_s0.2_t0.4.rds")
Wt <- c(cont_strong = 0, bin_strong = -0.0334, intermediate = 0.1154, near_null = 0.0496)
d <- do.call(rbind, lapply(names(cells), function(nm) { x <- readRDS(file.path(S, cells[[nm]]))$rows
  x$cell <- nm; x$wr <- (x$hi_pct - x$lo_pct) / (x$hi_wald - x$lo_wald)
  x$cov_wald <- x$lo_wald <= Wt[nm] & Wt[nm] <= x$hi_wald; x$cov_pct <- x$lo_pct <= Wt[nm] & Wt[nm] <= x$hi_pct; x }))
reg <- d[d$cell != "near_null", ]; nn <- d[d$cell == "near_null", ]
cat("Width ratio (percentile/Wald) quantiles, REGULAR cells pooled (n=600):\n")
print(round(quantile(reg$wr, c(.5, .9, .95, .99, .995, 1)), 3))
cat("per cell q95/q99/max:\n"); print(round(t(sapply(split(reg$wr, reg$cell), quantile, c(.95, .99, 1))), 3))
cat("\nNear-null width ratio quantiles:\n"); print(round(quantile(nn$wr, c(.1, .25, .5, .75, .9, 1)), 3))
cat("\nThreshold table: false-positive rate in regular cells | near-null flag rate | near-null flagged Wald cov | unflagged Wald cov | P(flag | oe_irregular) | P(flag | oe_regular)\n")
for (t in c(1.2, 1.25, 1.3, 1.5, 2, 3)) {
  f <- nn$wr >= t
  cat(sprintf("t=%.2f  FP=%.3f  NN_rate=%.3f  covW_flag=%.3f  covW_unfl=%.3f  P(flag|irreg)=%.3f  P(flag|reg)=%.3f\n",
      t, mean(reg$wr >= t), mean(f), mean(nn$cov_wald[f]), mean(nn$cov_wald[!f]),
      mean(f[nn$oe_regular == 0]), mean(f[nn$oe_regular == 1])))
}
cat("\nNear-null by A1: oe_regular=FALSE n=", sum(nn$oe_regular == 0), " Wald cov=", round(mean(nn$cov_wald[nn$oe_regular == 0]), 3),
    " pct cov=", round(mean(nn$cov_pct[nn$oe_regular == 0]), 3), " median Wald width/|W|=", round(median((nn$hi_wald - nn$lo_wald)[nn$oe_regular == 0]) / abs(Wt["near_null"]), 1), "\n")
cat("                 oe_regular=TRUE  n=", sum(nn$oe_regular == 1), " Wald cov=", round(mean(nn$cov_wald[nn$oe_regular == 1]), 3),
    " pct cov=", round(mean(nn$cov_pct[nn$oe_regular == 1]), 3), " median Wald width/|W|=", round(median((nn$hi_wald - nn$lo_wald)[nn$oe_regular == 1]) / abs(Wt["near_null"]), 1), "\n")
cat("\nA2 beyond A1 near-null: among oe_regular=TRUE draws, flagged(3x) n=", sum(nn$weak_id == 1 & nn$oe_regular == 1),
    " Wald cov flagged=", round(mean(nn$cov_wald[nn$weak_id == 1 & nn$oe_regular == 1]), 3),
    " unflagged=", round(mean(nn$cov_wald[nn$weak_id == 0 & nn$oe_regular == 1]), 3), "\n")
cat("\nPercentile vs Wald coverage by cell:\n")
print(round(t(sapply(split(d, d$cell), function(x) c(pct = mean(x$cov_pct), wald = mean(x$cov_wald),
  pct_width_over_wald = median(x$wr), boot_sd_cv = sd(x$se_boot) / mean(x$se_boot)))), 3))
cat("\nBinary: bootstrap SD outliers (se_boot/empSD > 2): n=", sum(d$se_boot[d$cell == "bin_strong"] / sd(d$W[d$cell == "bin_strong"]) > 2), "\n")
