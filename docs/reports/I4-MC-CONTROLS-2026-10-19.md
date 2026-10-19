# I4: Monte Carlo control comparison (2026-10-19)

Job 4339014: 180/180 tasks, 36,000 outer reps, 0 NA intervals, 0 failed resamples.
Gate and rationale: `tasks/plan-mc-bootstrap.md` (I4). Script: `inst/sim/i4_gate_analysis.R`;
per-cell table (4 candidates x 18 cells, coverage, paired diff, 90% CI, width ratio):
`inst/sim/results/calibration/i4_gate_results_2026-10-19.csv`.

## Verdict
No candidate passes every cell, so the prespecified rule gives no-go and retains (b).
Cells passed over each candidate's own applicable cells: a 2/18, b 16/18, c 16/18, d 11/12. (d is the reference in the 6 binary cells, so it is compared only in the 12 Gaussian cells.)

| Candidate | Failed cells |
|---|---|
| b CRN K=50k | s4/1000 (diff -0.007, CI [-0.0115,-0.0030]); s6/1000 (diff -0.011, CI [-0.0165,-0.0055]) |
| c variance correction | s4/1000 (-0.008); s6/250 (+0.0205) |
| d fresh K=200k (Gaussian cells only) | s6/1000 only (diff +0.0085, CI [+0.0055,+0.0120], width ratio 1.03) |
| a fresh K=5k | 16 of 18 cells (overcover up to 6 pp, widths up to 1.85x) |

## Truth source (changed after results were inspected)
The first pass used the 4-d.p. MC values in `truth.csv` (2e6 draws, SE ~3.5e-4). Null cells
are exactly 0.5, so the analysis was switched to closed form. This made b look worse
(s6/1000 diff -0.006 -> -0.011), i.e. against the candidate being retained.
Recomputing I2 coverage with closed-form truth reproduces its headline (s6/1000 0.944);
with MC truth the same rows give 0.9375. The manuscript pipeline (`code/analyze.R` in pmed)
reads `truth.csv`. The shift is cell-specific: closed-form truth lowered s6/1000 coverage by 0.65 pp but raised s7/1000 and s8/1000, so only s6/1000 is clearly affected. Paired
differences between methods are unaffected. Not edited; flagged for I7.

## Reproducibility
Runner and code are committed under `inst/sim/hopper/i4/`; sha256 match the I4 manifest.

## Mechanism for b's s6/1000 miss
Not noise. `docs/reports/CRN50K-RERUN-2026-10-14.md` (pmed) explains it: under CRN the local slope near
b = 0 is estimated from about 500 near-tie pairs at K = 50k, N = 1000, which randomly rescales each
interval, and the defect worsens with N. An independent-seed I6 is expected to reproduce it.
