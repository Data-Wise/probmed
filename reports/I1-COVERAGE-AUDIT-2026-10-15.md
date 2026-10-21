# Phase 7 I1: are failed replications silently dropped from the pmed coverage? (2026-10-15)

Question (plan `tasks/plan-mc-bootstrap.md`, Q4): pmed `code/analyze.R` computes
coverage with `na.rm = TRUE`, so a failed replication could vanish from the
denominator instead of counting as a miss. This audit recomputes coverage over
all 2,000 attempted replications per cell, with a missing or non-finite interval
counted as a miss. It is read-only on the pmed repo (`draft`).

## Result: no silent drops in any run

| Run | Folder | Cells (scenario × N × method) | Cells with ≠ 2000 rows | Non-finite intervals | max \|reported − all-attempted\| |
|---|---|---|---|---|---|
| fresh (original, K = 5k fresh draws) | `simulation_results_fresh/` | 70 | 0 | 0 | 0 |
| crn5k | `simulation_results_crn5k/` | 70 | 0 | 0 | 0 |
| crn50k (current) | `simulation_results/` | 70 | 0 | 0 | 0 |

No cell differs by more than 0.005, so none is flagged. Every reported pmed
coverage equals coverage over all attempted replications. **Q4 is answered:
the manuscript's coverage numbers are not inflated by dropped failures.**

**Positive control.** I planted a defect in the crn50k s1, N = 1000, MBCO cell:
40 intervals set to NA and 20 replications removed. Clean coverage is 0.947.
The planted reported coverage was 0.9464, against 0.918 over all attempted, a
gap of 0.028. The audit detects that kind of drop.

## Failure channels the output cannot show

The audit can only see rows that were written. Two channels sit inside the
simulation code:

1. **Failed bootstrap resamples** (`code/simulation_functions.R:118, 123`). A
   resample whose fit errors returns NA, and the percentile interval is taken
   with `quantile(..., na.rm = TRUE)`. The number of failed resamples per
   replication is not recorded, so an interval could rest on fewer than
   `n_boot = 1000` draws without trace. This cannot be audited from the saved
   results. Fix: write a `boot_failed` count per row in the next rerun (I6).
2. **Separated binary resamples are kept, not dropped.** `sep_warnings` counts
   resamples flagged for separation or non-convergence. They stay in the
   percentile distribution. crn50k, percentile bootstrap:

   | Cell | share of replications with any flagged resample | mean / max flagged resamples (of 1000) |
   |---|---|---|
   | s3, N = 50 | 0.865 | 135 / 994 |
   | s3, N = 100 | 0.890 | 68 / 981 |
   | s3, N = 250 | 0.381 | 4.0 / 319 |
   | s3, N = 500 | 0.023 | 0.05 / 18 |
   | s4, N = 50 | 0.047 | 0.11 / 45 |

   s3 (binary outcome) at N ≤ 100 builds most intervals partly from
   quasi-separated fits. The coverage there is honest, since nothing is
   dropped, but it is coverage of a bootstrap distribution that includes
   unstable fits. The manuscript should say so wherever it reports s3 at small
   N. MBCO rows exist only for Gaussian outcomes, so s3 and s4 are bootstrap-only.

Separately, a replication that errors entirely
(`code/run_simulation_chunk.R:49–51`) would be dropped as NULL. The row counts
show it never happened: 2,000 rows in every cell.

## Implications for Phase 7

- I2–I4 can use the stored coverages as computed. No re-scoring is needed.
- I6 (authoritative rerun) should record `boot_failed` per replication. It
  should also keep `sep_warnings` in the analysis table for s3 and s4.

## Script

Run from any directory. It reads `~/projects/research/pmed/simulation_results*/`.

```r
root <- "~/projects/research/pmed"
runs <- c(fresh = "simulation_results_fresh", crn5k = "simulation_results_crn5k",
          crn50k = "simulation_results")
REPS <- 2000L
out <- do.call(rbind, lapply(names(runs), function(r) {
  d <- read.csv(file.path(root, runs[[r]], "all_results.csv"))
  tr <- read.csv(file.path(root, runs[[r]], "truth.csv"))
  d <- merge(subset(d, estimator == "pmed"), tr[, c("scenario", "true_pmed")], by = "scenario")
  do.call(rbind, by(d, list(d$scenario, d$N, d$method), function(x) {
    ok <- is.finite(x$ci_lower) & is.finite(x$ci_upper)
    hit <- ok & x$ci_lower <= x$true_pmed & x$true_pmed <= x$ci_upper
    data.frame(run = r, scenario = x$scenario[1], N = x$N[1], method = x$method[1],
               rows = nrow(x), unique_reps = length(unique(x$rep)), finite_ci = sum(ok),
               cov_reported = mean(x$ci_lower <= x$true_pmed & x$true_pmed <= x$ci_upper, na.rm = TRUE),
               cov_all = sum(hit) / REPS)
  }))
}))
out$diff <- out$cov_reported - out$cov_all
subset(out, abs(diff) > 0.005 | finite_ci < REPS)   # 0 rows on 2026-10-15
```
