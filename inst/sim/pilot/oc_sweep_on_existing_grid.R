## Threshold sweep run on the EXISTING 16,000-rep coverage grid, using
## se_ratio = seW_bt/seW_an as a stand-in for the percentile/Wald CI-width ratio.
##
## READ THIS BEFORE CITING ANY NUMBER THIS PRODUCES. The stand-in does NOT track
## the shipped weak_id_ratio (r ~ 0.41; see se_ratio_proxy_check.R), so the
## "optimum near 3" this prints is NOT a validation of the shipped
## weak_id_ratio_threshold = 3. It is retained because (a) it documents the
## attempt, so nobody re-derives the shortcut and mistakes it for evidence, and
## (b) its metric-INDEPENDENT output is genuinely informative: overall Wald
## coverage ~0.88 across all 16k reps, i.e. sub-nominal regardless of any flag.
##
## Requires ~/gauge_boot/gauge_boot_raw.rds on hopper (not committed -- 660KB of
## raw reps; the collated per-cell CSV is in ../results/).
## Run on hopper:  Rscript oc_sweep_on_existing_grid.R
d <- readRDS(path.expand("~/gauge_boot/gauge_boot_raw.rds"))
d$se_ratio <- d$seW_bt / d$seW_an
cat(sprintf("reps: %d\n", nrow(d)))
cat(sprintf("se_ratio: median %.2f  IQR [%.2f, %.2f]  max %.1f\n",
            median(d$se_ratio, na.rm = TRUE),
            quantile(d$se_ratio, .25, na.rm = TRUE),
            quantile(d$se_ratio, .75, na.rm = TRUE),
            max(d$se_ratio, na.rm = TRUE)))
cat("\nDOES A RATIO-THRESHOLD FLAG PREDICT BAD WALD COVERAGE? (PROXY -- see header)\n")
cat(sprintf("%8s %9s %11s %13s %11s\n",
            "thresh", "pct_flag", "covW_flag", "covW_unflag", "separation"))
for (t in seq(1.5, 6, by = 0.5)) {
  fl <- d$se_ratio >= t
  fl[is.na(fl)] <- FALSE
  if (sum(fl) < 50 || sum(!fl) < 50) next
  cf <- mean(d$covW_an[fl], na.rm = TRUE); cu <- mean(d$covW_an[!fl], na.rm = TRUE)
  cat(sprintf("%8.2f %9.3f %11.3f %13.3f %11.3f\n", t, mean(fl), cf, cu, cu - cf))
}
cat(sprintf("\nMETRIC-INDEPENDENT (this one IS citable): overall Wald coverage %.3f\n",
            mean(d$covW_an, na.rm = TRUE)))
