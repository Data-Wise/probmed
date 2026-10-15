# Calibration results: single-mediator repair (T16–T17), 2026-10-14

Prespecified design: pmed `docs/SPEC-probmed-repair-2026-10-09.md`; harness
`inst/sim/calibration/` at `af243ab`. 36 cells × 2000 replications, base seed
20261014, L'Ecuyer-CMRG stream per replication. Run locally (14 cores, ~52
core-h) because hopper has medfit 0.3.1 and no lpSolve.

Results: `inst/sim/results/calibration/` (`calibration.csv.gz` per replication,
`calibration_summary.csv` per cell, `calibration_stability.csv` 1× vs 4× controls).

Decision rule (frozen): an exact 95% binomial CI for coverage over **all**
attempted replications, where an incomplete interval counts as a miss. PASS if
the CI lies inside [0.93, 0.97]; FAIL if it lies wholly outside; INCONCLUSIVE
otherwise.

**Totals: 17 PASS · 15 FAIL · 4 INCONCLUSIVE.**

## MBCO (`method = "mbco"`), 21 cells

| Pair (a, b) | N=50 | N=250 | N=1000 |
|---|---|---|---|
| (0, 0) double null | FAIL over 0.997 | FAIL over 0.997 | FAIL over 0.998 |
| (0, 0.4) | PASS 0.958 | PASS 0.953 | PASS 0.946 |
| (0.5, 0) | INCONCL 0.973 | PASS 0.950 | PASS 0.947 |
| (0.05, 0.05) near null | FAIL over 0.994 | FAIL over 0.994 | INCONCL 0.976 |
| (0.5, 0.4) | PASS 0.949 | PASS 0.944 | PASS 0.954 |
| (−0.5, 0.4) | PASS 0.952 | PASS 0.949 | PASS 0.949 |
| (2.5, 1.5), Vy=2 | PASS 0.947 | PASS 0.945 | PASS 0.946 |

All 42,000 MBCO intervals are complete, so conditional coverage equals overall
coverage. Every FAIL is **overcoverage**, and all of them sit at the composite
null a=b=0 or next to it. This is the expected price of the exact-null union
{a=0}∪{b=0}: the LR statistic at the double null is not χ²₁, so the interval is
conservative there. A single-null cell (one path nonzero) and every non-null
cell is nominal.

**Supported claim:** nominal coverage at single-null and non-null
configurations. The interval is conservative (≈0.997) at and near the double null.

## Model bootstrap (`method = "model_bootstrap"`), 15 cells

| Config | N=50 | N=250 | N=1000 |
|---|---|---|---|
| 1 Gaussian M, Gaussian Y | PASS 0.950 | PASS 0.961 | INCONCL 0.977 |
| 2 Gaussian M, binary Y (5%) | FAIL under 0.000 | FAIL under 0.227 (cond 0.981) | FAIL over 0.985 |
| 3 Gaussian M, binary Y (50%) | FAIL under 0.910 (cond 0.964) | PASS 0.958 | FAIL over 0.983 |
| 4 binary M, Gaussian Y | INCONCL 0.937 (cond 0.994) | FAIL over 0.987 | FAIL over 0.999 |
| 5 binary M, binary Y, + Z | FAIL under 0.842 (cond 0.998) | FAIL over 0.993 | FAIL over 0.999 |

There are three distinct failure mechanisms:

1. **Fit-failure strata (rare outcome, small N).** Config 2 at N=50 never
   yields a complete interval: the original fit is separated, or the
   bootstrap refits fail past the 99% completion rule. At N=250, 77% of
   replications are incomplete. Coverage *among complete intervals* is 0.981.
   These cells measure how often the method declines to report, not
   miscalibration of the intervals it does report.
2. **Small-N incompleteness.** In configs 3–5 at N=50, undercoverage comes
   entirely from incomplete intervals, which are counted as misses. Conditional
   coverage is 0.964–0.998.
3. **Large-N overcoverage, consistent with Monte Carlo inflation.** Each
   bootstrap draw evaluates P_med by simulation with `n_sim = 1e4`. That
   simulation noise (sd ≈ 0.004–0.006, measured below) is added to the
   bootstrap spread. At N=1000 it is comparable to the sampling SE of the
   estimator, so the percentile interval widens.

### Stability subset (100 reps per cell, 4× `n_boot` and 4× `n_sim`)

| Cell | cov 1× | cov 4× | width ratio 4×/1× |
|---|---|---|---|
| boot-1-N1000 | 1.00 | 0.97 | 0.89 |
| boot-2-N1000 | 1.00 | 0.98 | 0.88 |
| boot-3-N1000 | 0.98 | 0.96 | 0.87 |
| boot-4-N1000 | 1.00 | 0.98 | 0.75 |
| boot-5-N1000 | 1.00 | 1.00 | 0.73 |
| boot-4-N250 | 0.99 | 0.94 | 0.91 |
| boot-5-N250 | 0.99 | 0.97 | 0.90 |
| N=50 cells | — | — | 0.98–1.01 |

The width ratio is about 1 at N=50 and drops monotonically with N, to 0.73 for
config 5 at N=1000. Coverage moves toward nominal at the same time. This is the
signature of an additive MC-variance term that does not shrink with N. The 4×
controls do not fully remove it at N=1000, and 100 replications cannot
separate 0.95 from 0.97. **The hypothesis is supported, not confirmed.** The
next step would be `n_sim` scaled with N (for example, MC sd below 0.25 ×
the estimated SE) and a rerun of the N ≥ 250 boot cells.

**Supported claim:** nominal or conservative coverage for complete intervals at
N ≥ 250 for configs 1, 3, 4 and 5. Overcoverage at large N is driven by the
simulation noise in `n_sim`. Rare binary outcomes (5%) at N ≤ 250 are a
fit-failure regime: the method frequently returns `bootstrap_incomplete`
rather than an interval, and users must be told so.

## What this does not support

- A blanket "95% coverage" claim for `model_bootstrap`. Only configs 1 and 3 at
  N=250 PASS outright.
- Any MBCO claim at the double null beyond "conservative".
- Generalization outside this design: single mediator, binary X, with or without
  one binary covariate.

## Addendum (2026-10-14): T20 positive control, Monte Carlo cause confirmed

`inst/sim/calibration/calib_t20_exact.R` reran boot-1-N1000 (Gaussian M and Y)
with the main run's streams, so all 2000 datasets are identical to the main run.
The only change was the exact P_med inside every bootstrap draw in place of the
`n_sim = 1e4` simulation. Before the swap, a known-answer check confirmed the
exact formula matches the simulation (0.5522 vs 0.5522 at n_sim = 4e6).
Results: `inst/sim/results/calibration/t20_exact_boot1_N1000.csv`.

| | Coverage [exact 95% CI] | Median width | Decision |
|---|---|---|---|
| Main run (MC in each draw) | 0.977 [0.969, 0.983] | 0.0354 | INCONCLUSIVE |
| Exact P_med in each draw | 0.945 [0.934, 0.955] | 0.0296 | PASS |

- Removing the simulation noise narrows the interval by 16% and brings
  coverage to nominal.
- Every disagreement goes one way: 64 replications covered only with the
  simulation and none only with the exact version.
- The simulation error in the point estimate has sd 0.0049, against the
  estimator's own sampling sd of 0.0076.

The overcoverage is therefore Monte Carlo inflation, not a defect of the
response bootstrap. **Implication for T21:** common random numbers remove the
noise from the percentile spread. The point estimate still carries a fixed
simulation error of about 0.005 at `n_sim = 1e4`, about two-thirds of the
sampling sd at N = 1000. T21 should therefore also use a larger `n_sim` for the
point estimate alone, which costs one simulation, not `n_boot` of them.

## Addendum (2026-10-14): T22, prespecified rerun of the bootstrap cells after T21

Design: `inst/sim/calibration/calib_manifest_t22.dcf`, committed before the full
array ran. The run covers the same 15 cells and 2000 replications per cell,
base seed 20261015, with the frozen decision rule unchanged. It used probmed
`41d71de`: common random numbers, draws re-centered on a control variate, and a
point estimate with `100 * n_sim` pairs. It ran on hopper (array 4337068, 300
tasks, all completed) from a fresh library with medfit 0.3.2 and lpSolve 5.6.23.

One pilot task at the earlier commit `4bba493` (uncentered shared uniforms)
covered 0.79 in boot-5-N1000. The cause was the common Monte Carlo error
shifting the whole interval, and it was fixed before the array ran. That pilot
output is set aside and not collated. Results:
`inst/sim/results/calibration/t22/`.

**Totals: 9 PASS · 4 FAIL · 2 INCONCLUSIVE** (the original run had 2 · 10 · 3 on
these cells).

| Config | N=50 | N=250 | N=1000 |
|---|---|---|---|
| 1 Gaussian M, Gaussian Y | PASS 0.946 | PASS 0.946 | INCONCL 0.934 [0.923, 0.945] |
| 2 Gaussian M, binary Y (5%) | FAIL 0.000 (no complete intervals) | FAIL 0.242 (cond 0.955) | PASS 0.943 |
| 3 Gaussian M, binary Y (50%) | FAIL 0.904 (cond 0.951) | PASS 0.947 | PASS 0.946 |
| 4 binary M, Gaussian Y | INCONCL 0.931 (cond 0.990) | PASS 0.952 | PASS 0.942 |
| 5 binary M, binary Y, + Z | FAIL 0.847 (cond 0.999) | PASS 0.949 | PASS 0.956 |

- **The large-N overcoverage is gone.** All seven cells that failed for
  overcoverage now PASS. Widths shrank by up to 37%, for example boot-5-N1000
  from 0.0168 to 0.0106.
- **Every remaining FAIL is a fit-failure or incompleteness stratum.** Coverage
  among complete intervals is 0.951–0.999 in each of them. Config 2 at N ≤ 250
  is too sparse to fit reliably. The N=50 cells in configs 3–5 lose intervals
  to the 99% completion rule.
- **boot-1-N1000 is INCONCLUSIVE** under the frozen rule: 0.9345, exact 95% CI
  [0.9228, 0.9449], which excludes 0.95. T20's exact per-draw evaluation of the
  same design gave 0.945. Remaining point-estimate Monte Carlo error or a
  fresh seed may explain the gap, but that is an untested inference, not a
  result. The cell needs a dedicated rerun before any nominal claim (see the
  Phase 7 plan, `tasks/plan-mc-bootstrap.md`).

**Supported claim, replacing the T17 claim for `model_bootstrap`:** nominal
coverage (PASS) at N ≥ 250 in configs 2–5, except a rare (5%) binary outcome
at N = 250, and in config 1 at N = 250. **Config 1 (Gaussian M and Y) at
N = 1000 is INCONCLUSIVE and leans low (0.9345).** At N = 50, and for rare outcomes at
N ≤ 250, the method often returns `bootstrap_incomplete`. The intervals it
does report cover at nominal or conservative rates.

*Corrected 2026-10-14 after the Codex adversarial review: the first version of
this addendum called boot-1-N1000 nominal.*

## Addendum (2026-10-14): MBCO at the composite null (Phase 8, T2/T4)

Harness: `inst/sim/calibration/calib_mbco_null.R` at `e718a26`, hopper arrays
4338241 (pilot) and 4338242. A √N-local lattice (a, b) = (h_a, h_b)/√N with
h ∈ {0, 0.5, 1, 2, 4}², N ∈ {100, 1000}, with or without a covariate:
100 points × 2,000 replications. Self-check: interval coverage equals
acceptance of the truth in 5,000 of 5,000 full-inversion replications, with 0
failed profiles in 200,000. Results: `inst/sim/results/calibration/mbco_null/`.

**χ²₁ (shipped) never undercovers.** Its minimum lattice coverage is 0.9415
(N = 100) and 0.9440 (N = 1000), within Monte Carlo error of 0.95. It is
conservative across a band, not only at the origin: coverage is 0.97–0.999
wherever both h ≤ 2 (both paths within about 2 SE of zero), and 0.997–0.999
at a = b = 0. The new `near_double_null` diagnostic (both branch LR below the
cutoff, i.e. both |z| < 1.96) flags this band.

**No simple calibration is acceptable.** The prespecified acceptance was
minimum coverage ≥ 0.94 and coverage at the origin ≤ 0.97.

| Rule | min coverage N=100 / N=1000 | origin N=100 / N=1000 |
|---|---|---|
| χ²₁ | 0.9415 / 0.9440 | 0.997 / 0.9985 |
| switch to qchisq(1 − √α) when both branches < χ² cutoff | 0.863 / 0.875 | 0.966 / 0.965 |
| same, threshold log N | 0.854 / 0.836 | 0.958 / 0.952 |
| smooth interpolation on max branch LR | 0.908 / 0.917 | 0.997 / 0.9985 |

Lowering the cutoff fixes the origin but undercovers at nearby non-null points
such as h = (4, 1) and (2, 2). There P_med ≠ 0.5 and LR(truth) is no longer the
minimum of the branch statistics. This is the non-uniformity that the
composite-null literature addresses (Miles & Chambaz; van Garderen & van
Giersbergen; He, Song & Xu). A valid calibrated interval would need one of
those constructions, inverted over all p*, not only the p* = 0.5 test.

**Decision for now:** keep the χ²₁ interval. Label it conservative near the
double null (PR #49), and do not ship a calibrated cutoff.
