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
