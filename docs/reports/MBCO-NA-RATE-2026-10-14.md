# MBCO interval status rate after the repair (task T3)

Date: 2026-10-14. Script: `inst/sim/mbco_na_rate.R` (200 seeded datasets per cell,
16 cells: effect {a = 0, b = 0, weak 0.15/0.15, strong 0.5/0.4} × n {100, 500} ×
{no covariate, one covariate}; Gaussian M and Y, continuous X). Result:
`inst/sim/results/mbco_na_rate.csv`. R 4.6.1, medfit 0.7.0.

## Before the fix (batch 2 as received)

P_med intervals came back `failed` in 1–2.5% of datasets in the weak cells and in
the strong / n = 100 / covariate cell (1.5%), above the 1% threshold for a
well-conditioned cell. Every failure had reason `root_verification`.

Cause: near the null, each sign basin of the signed P_med profile has two branches,
a ≈ 0 (b near b̂) and b ≈ 0 (a near â). Both optimizer starts sat on the a ≈ 0
branch. For p just above 0.5 they stayed there although the b ≈ 0 branch has the
higher likelihood (seed 1200018: −296.2 vs −292.4), so the profile jumped by ~6.7
between p = 0.50066 and 0.50067. The root check correctly refused the jump; the
interval was reported missing rather than wrong.

## Fix

`.mbco_fit_pmed()` adds a third start aimed at the b ≈ 0 branch: holding a at |â|
in the constraint gives r = (â·δ/q*)²/2 − V̂m. Regression test
"near-null P_med profiles reach the b ~ 0 branch of each sign basin" in
`tests/testthat/test-repair-profile.R` fails without the fix (3 failures, 1 error)
and passes with it, agreeing with the independent nloptr oracle within 1e-6 per
observation and placing the lower endpoint on the χ²₁ critical value within 1e-3.

## After the fix

| | P_med `failed` | IE `failed` |
|---|---|---|
| Max over 16 cells | 0 / 200 | 1 / 200 (0.5%) |
| Cells with any failure | 0 | 5 |

No `truncated` or `disconnected` statuses occurred. Median wall time per call
(P_med + IE intervals) was 2.0–2.6 s, 95th percentile ≤ 3.0 s.

Thresholds (plan, 2026-10-14): ≤ 1% well-conditioned, ≤ 5% near the null. Met.

These are status rates only. They say nothing about coverage; calibration is a
separate deliverable (tasks T15–T17).
