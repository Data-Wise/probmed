## Does se_ratio = seW_bt/seW_an track the SHIPPED weak_id_ratio (a percentile/Wald
## CI-WIDTH ratio)?  (Cited in ../hopper/README.md and .STATUS item 5.)
##
## WHY THIS EXISTS: the completed 16k-rep coverage grid stores no CI widths, so it
## is tempting to run the #11 threshold sweep on se_ratio as a stand-in -- and that
## sweep looks encouraging (separation peaks near 3, apparently confirming the
## shipped default). This script is the check that kills that shortcut: the proxy
## correlates with the shipped metric at only r ~ 0.41 (Spearman ~0.53), on a
## different scale (medians ~2.4 vs ~1.7), agreeing on the ">= 3" call ~72% of the
## time. A threshold calibrated on se_ratio therefore says little about the same
## number on the width scale.
##
## Caveat on this script's own evidence: n = 64. The Fisher-z 95% CI for r = 0.41
## is roughly [0.18, 0.60] -- enough to rule OUT the strong correlation a transfer
## would need, even at the upper bound, but not a precise estimate of r.
##
## Run:  Rscript inst/sim/pilot/se_ratio_proxy_check.R    (64 fit-pairs; ~15 min)
suppressMessages(devtools::load_all(quiet = TRUE))
expit <- function(x) 1 / (1 + exp(-x))

gen <- function(n, s, binY, seed) {
  set.seed(seed)
  C <- rnorm(n); A <- rbinom(n, 1, expit(-0.2 + 0.8 * C))
  M <- s * 0.6 * A + 0.4 * C + rnorm(n)
  lin <- s * 0.5 * A + 0.7 * M + s * 0.4 * A * M + 0.3 * C
  Y <- if (binY) rbinom(n, 1, expit(lin)) else lin + rnorm(n)
  data.frame(A, M, Y, C)
}

## span the identification range so r is not measured within one regime
grid <- expand.grid(s = c(0.05, 0.15, 0.30, 1.0), n = c(800, 3000))
res <- list(); k <- 0L
for (g in seq_len(nrow(grid))) {
  for (i in 1:8) {
    k <- k + 1L
    d  <- gen(grid$n[g], grid$s[g], FALSE, 5000 + k)
    fa <- ward_residual(d, seed = 5000 + k)
    fb <- suppressWarnings(
      ward_residual(d, seed = 5000 + k, se_method = "bootstrap", B = 200L))
    res[[k]] <- data.frame(s = grid$s[g], n = grid$n[g],
      width_ratio = fb@weak_id_ratio,     # SHIPPED metric
      se_ratio    = fb@W_se / fa@W_se,    # the stand-in
      oe_snr      = fb@oe_snr)
  }
}
r <- do.call(rbind, res)
r <- r[is.finite(r$width_ratio) & is.finite(r$se_ratio), ]
cat(sprintf("n reps: %d\n", nrow(r)))
cat(sprintf("Pearson  r = %.3f\n", cor(r$width_ratio, r$se_ratio)))
cat(sprintf("Spearman r = %.3f  (rank agreement -- what matters for a threshold)\n",
            cor(r$width_ratio, r$se_ratio, method = "spearman")))
cat(sprintf("\nmedian width_ratio %.2f  vs  median se_ratio %.2f\n",
            median(r$width_ratio), median(r$se_ratio)))
cat("\n--- do they agree on WHICH reps exceed 3? ---\n")
tab <- table(width_ge3 = r$width_ratio >= 3, se_ge3 = r$se_ratio >= 3)
print(tab)
cat(sprintf("\nagreement: %.1f%%\n", 100 * sum(diag(tab)) / sum(tab)))
