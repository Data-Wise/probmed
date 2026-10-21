# Observed-data Wasserstein validation: merge blocked

Source: d476a0f10ec76993dda05e3c7cfd2eb890149751. Design and frozen
diagnostics: [repair plan](PLAN-wasserstein-observed-inference-2026-10-20.md).
The repair has 35 focused assertions and 42 existing Wasserstein assertions
passing. These are not a full-package local-CI result.

## Repeated sampling

Each cell has 2,000 independent datasets, zero failed fits. Bias MCSEs are
0.001266, 0.001542, 0.000890, and 0.000994 in table order. The absolute
bias/SD threshold is <0.10, SE/SD must be in [0.90,1.10], and coverage in
[0.92,0.98]. No threshold has been changed after observing results.

| Outcome | n | Bias | Bias/SD | Mean SE/SD | Coverage | Verdict |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Gaussian | 500 | -0.001286 | 0.0227 | 1.00298 | 0.9485 | Pass |
| Binary | 500 | 0.009851 | 0.1428 | 0.98839 | 0.9630 | **Fails bias diagnostic** |
| Gaussian | 1,000 | -0.001241 | 0.0312 | 1.00879 | 0.9490 | Pass |
| Binary | 1,000 | 0.004159 | 0.0935 | 0.97302 | 0.9635 | Pass |

The observed-data uncertainty repair addresses the original simulation-count
defect, but acceptable SE calibration and coverage do not erase a failed bias
diagnostic. On 2026-10-20 the author explicitly retained the merge hold and
requested investigation of bias correction. Main protection and hosted CI
must remain unchanged while the scientific issue is unresolved.

## Reproducibility and limits

Run `Rscript scripts/validate-wasserstein-observed.R OUTPUT_DIRECTORY` from the
package root. Saved [replicates](data/wasserstein-observed-2026-10-20/replicates.csv),
[summary](data/wasserstein-observed-2026-10-20/summary.csv), and
[provenance](data/wasserstein-observed-2026-10-20/provenance.txt) preserve all
8,000 outcomes, seeds, versions, and elapsed time (274.08 seconds).
SHA-256 of replicates.csv:
`91f65e70b837f6e4fcabbf1f76f5722687a791b1723685426c6b1e4f7d9212a6`.

Finite-difference weighted-refit tests support the influence calculation;
omitting empirical covariate uncertainty causes a targeted test failure.
Independent Gaussian and binary integration checks support selected numerical
fixtures. These are not proofs of uniform accuracy, efficiency, boundary
coverage, or validity under misspecified models. Fixed-order quadrature error
must remain negligible relative to sampling error for the approximation used.
The legacy `eif` name denotes model-based observed-data inference, not a
nonparametric-efficiency claim.

The [binary bias investigation](INVESTIGATE-wasserstein-binary-bias-2026-10-20.md)
records the completed decomposition and the exploratory correction design.
The Wasserstein ratio is not the original tie-adjusted P_med estimand. No
P_med/I6/I7, release, or submission gate is changed by this investigation.
