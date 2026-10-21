# Findings: the bootstrap arm and the two gates, re-measured after the weight fix

**Date:** 2026-08-22 (evening), after PR #34 (`a73efc8`) corrected the corner-EIF
inverse-probability weights.
**Script:** `inst/sim/boot_gates_check.R` — `ward_residual(se_method = "bootstrap",
B = 200, reps = 1)` on 200 datasets per cell, n = 800, against the exact truth of
the `phi_decomposition.R` DGPs, oracle-nuisance se as benchmark. Calibration
analysis `inst/sim/boot_gates_calibrate.R` over the saved per-dataset rows;
per-cell summaries in `inst/sim/results/boot_gates_postfix.csv`.
**Question:** PR #34 left the bootstrap arm and the `weak_id` (A2) / `oe_regular`
(A1) gates unmeasured with correct weights. What do they do now, and can A2's
3x threshold be recalibrated?

## Cells

| cell | `s` | `tau` | Y | `oe_snr` median | percentile cov | Wald cov | oracle cov |
|---|---|---|---|---|---|---|---|
| continuous, strong ID | 1 | 0 | cont | 9.5 | 0.905 | 0.935 | 0.940 |
| binary, strong ID | 1 | 0 | bin | 4.9 | 0.930 | 0.970 | 0.965 |
| intermediate | 0.5 | 0.4 | cont | 5.1 | 0.920 | 0.910 | 0.965 |
| near-null | 0.2 | 0.4 | cont | 1.9 | 0.990 | 0.945 | 0.990 |

Monte-Carlo half-width on a coverage at 200 datasets: about +/- 0.03.

## 1. The bootstrap (percentile) arm

| cell | boot SD / empSD (median) | CV of boot SD | percentile width / calibrated | percentile / Wald width (median, q90) |
|---|---|---|---|---|
| continuous, strong | 0.92 | 0.19 | 0.92 | 0.99, 1.11 |
| binary, strong | 1.03 | **2.18** | 1.01 | 1.05, 1.21 |
| intermediate | 1.00 | 0.52 | 0.99 | 1.02, 1.18 |
| near-null | 2.92 | 3.74 | 1.04 | 1.74, 3.88 |

- In the regular cells the percentile interval is **as wide as the Wald interval**
  (ratio median 0.99-1.05) and covers 0.905-0.930: at or slightly below nominal,
  never above. The pre-fix "~1.00 in every cell" is gone; so is the ~2x width gap
  the A2 threshold was set against.
- For binary Y the refit bootstrap is unstable: 11/200 datasets have a bootstrap SD
  more than 2x the empirical SD (CV 2.2), where the analytic se's CV is 0.33. A
  resample that lands a near-zero `OE_hat` dominates.
- Near the null the bootstrap SD explodes (2.9x, CV 3.7) exactly as the oracle's
  does: the ratio is non-regular and no construction fixes it. Both arms
  over-cover there (0.99 / 0.945).

**Reading:** the bootstrap arm no longer has a job. It is not more conservative,
not better calibrated, costs B x reps refits, and is unstable for binary Y. Keep
it as an option for users who want a refit bootstrap; stop presenting it as the
remedy for anything.

## 2. A2 — `weak_id` (percentile / Wald width ratio >= 3)

Width-ratio distribution, regular cells pooled (n = 600): median 1.02, q90 1.17,
q95 1.22, q99 1.49, max 2.79 (one binary outlier). Near-null: median 1.74, q90
3.88, max 9.2.

Threshold sweep (FP = flag rate in the regular cells; near-null columns from the
200 near-null datasets):

| t | FP | near-null flag rate | Wald cov, flagged | Wald cov, unflagged | P(flag \| A1 irregular) | P(flag \| A1 regular) |
|---|---|---|---|---|---|---|
| 1.20 | 0.068 | 0.83 | 0.940 | 0.971 | 0.88 | 0.78 |
| 1.25 | 0.040 | 0.80 | 0.938 | 0.975 | 0.88 | 0.71 |
| 1.30 | 0.030 | 0.75 | 0.953 | 0.922 | 0.86 | 0.62 |
| 1.50 | 0.010 | 0.62 | 0.951 | 0.935 | 0.82 | 0.40 |
| 2.00 | 0.002 | 0.41 | 0.951 | 0.941 | 0.65 | 0.15 |
| 3.00 | 0.000 | 0.19 | 0.947 | 0.944 | 0.35 | 0.02 |

- **At no threshold do flagged draws under-cover.** Flagged and unflagged Wald
  coverage are within Monte-Carlo noise of each other at every t. The ratio is
  not a coverage predictor; it cannot be calibrated to one.
- What it does track is **bootstrap instability**, which is what A1 already
  detects: P(flag | `oe_regular = FALSE`) is 2-17x P(flag | TRUE) at t >= 1.5.
  Beyond A1, at the shipped 3x threshold, A2 flags **2 of the 98** near-null draws
  that pass A1.
- The one under-covering subset in the whole study is *near-null draws that pass
  A1* (Wald 0.898, n = 98; their intervals are a median 4.1x |W| wide). The width
  ratio does not separate them either (above: flagged vs unflagged coverage).

**Reading:** A2 cannot be recalibrated because there is no signal to calibrate it
to. Its pre-fix validation ("separation explained by Wald width", "narrowness
proxy") described the weight bug's 30% se shortfall, which no longer exists.
Recommended disposition: stop *warning* on `weak_id`; keep the field and ratio
computed (cheap, and a user running the bootstrap may want to see its
instability), documented as "bootstrap-vs-Wald width discrepancy; not a coverage
diagnostic; redundant with `oe_regular`". Leave the threshold at 3 rather than
moving it — a new number would imply a calibration that does not exist.

## 3. A1 — `oe_regular` (`oe_snr >= 2`)

Fires only in the near-null cell (51%; `oe_snr` median 1.9; Fieller bounded 49%).
Flagged draws **over**-cover (Wald 0.99, percentile 1.00; median Wald width
8.7x |W|); unflagged near-null draws cover 0.898 (width 4.1x |W|). Same reading as
the pre-fix grid: `oe_regular = FALSE` means "uninformative interval", not
"missed". The non-regularity is the estimand's, so the bug did not touch it.

**Reading:** A1 stands as is. The residual under-coverage among near-null draws
that pass it (0.90) is the local-to-zero regime no sample gate separates; the
Fieller set is the honest interval there, and the docs already say so.

## What this changes in the package

1. `weak_id`: drop the `warning()`; rewrite the property docs (A2 is a
   bootstrap-instability discrepancy, redundant with A1, not a coverage
   diagnostic; the 48k-grid "validation" paragraphs describe the pre-fix estimator).
2. `se_method` docs: the bootstrap arm is an option, not a remedy; percentile
   coverage 0.905-0.930 in regular cells; unstable for binary Y.
3. `oe_regular`: docs unchanged in substance; replace "measured pre-fix, expected
   to stand" with these numbers.
4. Promote `boot_gates_check.R` + `calibrate.R` to `inst/sim/`, results CSV to
   `inst/sim/results/`.

Code change: needs a worktree.

## Not measured

n = 3000; binary near-null; `reps > 1` under the bootstrap; the 48k gate grid
itself (hopper). None of these is expected to change the readings above: the
regular-cell width ratios are already tight (q99 1.49) and the near-null behavior
is oracle-limited.
