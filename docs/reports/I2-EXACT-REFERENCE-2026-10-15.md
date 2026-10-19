# I2: exact-bootstrap Gaussian reference grid (2026-10-15)

## Result

The exact percentile-bootstrap reference was computed for s1, s2, s6, s7, and s8 at N = 50, 100, 250, 500, and 1,000, with 2,000 replications and 1,000 resamples per replication. It uses the same datasets and bootstrap-resample stream as the pmed `crn50k` run. Across all 50,000 replications, the exact point estimate matched the saved MBCO point estimate within `4.0e-15`; this validates dataset replay. Bootstrap resamples were replayed by preserving the historical task seed, L'Ecuyer stream, 8-worker `mclapply` chunking, and random draws consumed before bootstrap sampling. There were no failed bootstrap resamples.

For each Gaussian resample, the exact fitted-model score is

`Phi(a_hat * b_hat / sqrt(2 * (b_hat^2 * sigma_M_hat^2 + sigma_Y_hat^2)))`,

from the normal difference between the two nested outcomes under independent outcome errors. The script replaces only simulated evaluation of each resample's score with the closed form.

| Scenario | N=50 | N=100 | N=250 | N=500 | N=1000 |
|---|---:|---:|---:|---:|---:|
| s1 | .9390 | .9545 | .9450 | .9540 | .9445 |
| s2 | .9375 | .9545 | .9555 | .9405 | .9525 |
| s6 (b=0) | .9830 | .9670 | .9405 | .9455 | .9440 |
| s7 (a=0) | .9590 | .9550 | .9565 | .9415 | .9475 |
| s8 | .9340 | .9420 | .9465 | .9410 | .9470 |

Each cell has 2,000 attempted replications; exact-binomial 95% intervals are in `inst/sim/results/calibration/i2_exact_summary_2026-10-15.csv`. The prespecified positive control passes: s6, N = 1,000 covers in 1,888/2,000 replications (0.944; exact-binomial 95% CI 0.9330–0.9537), which contains 0.951. It is consistent with 0.951; it does not equal it.

At the b=0 boundary, exact-bootstrap coverage is high at N = 50 and 100 (.983 and .967), then .941–.946 at N ≥ 250. The other composite-null branch (s7, a=0) ranges .942–.959 across the grid. These results support near-nominal coverage from N = 250 onward for the two tested Gaussian null mechanisms, but show small-sample overcoverage for s6; they do not establish coverage for other models or controls.

Exact mean widths match the saved CRN widths closely: exact/CRN width ratios range 0.9988–1.0003 across 25 cells. Exact widths are 0.9658–0.9981 times the saved MBCO widths (MBCO is about 0.2%–3.5% wider). MBCO coverage ranges 0.9400–0.9815; the largest absolute paired exact-versus-MBCO coverage difference is 0.0125 (s8, N = 50: 0.9340 versus 0.9465). This is a descriptive comparison, not a method-selection rule.

## Reproduction and provenance

The exact run was launched from the read-only pmed repository:

```sh
cd /Users/davood/projects/research/pmed
PMED_SEED_BASE=20261016 I2_CORES=8 Rscript --vanilla /private/tmp/probmed_i2_grid.R
```

The pmed run used seed base `20261016` (`docs/reports/CRN50K-RERUN-2026-10-14.md`) and source commit `1c0f46b`. At execution, pmed `HEAD` was `c6c3385`; the four source files used by this computation (`code/config.R`, `code/simulation_functions.R`, `code/mbco_pmed.R`, and `code/run_simulation_chunk.R`) were unchanged from `1c0f46b`. The pilot first used the wrong seed base `20260605` and differed from historical MBCO point estimates by as much as 0.166; using the documented seed base and original worker topology reduced the pilot discrepancy below `5e-15`. The complete run's maximum discrepancy was `4.0e-15`.

Files: `inst/sim/results/calibration/i2_exact_reference_2026-10-15.csv` (50,000 row-level intervals), `inst/sim/results/calibration/i2_exact_summary_2026-10-15.csv` (25 cell summaries), and `inst/sim/results/calibration/i2_exact_manifest_2026-10-15.csv` (input/output hashes and run settings). The run command was `PMED_SEED_BASE=20261016 I2_CORES=8 Rscript --vanilla /private/tmp/probmed_i2_grid.R` from the pmed root. The one-off runner stayed in `/private/tmp` because the active `dev` branch blocks new code files. The pmed repository was read-only; its pre-existing untracked files were left untouched.

## Interpretation and next gate

I2 answers Q1 within the five Gaussian scenarios and tested sample sizes. Exact percentile-bootstrap coverage is close to 95% from N = 250 onward in both null branches, while s6 shows small-sample overcoverage. I4 is next only after prespecifying and agreeing on its paired coverage-equivalence margin and rationale; do not infer equivalence from overlapping intervals.
