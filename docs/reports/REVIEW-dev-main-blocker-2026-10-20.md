# Dev-to-main review: Wasserstein inference blocker

Date: 2026-10-20 (America/Denver). Reviewed head: `5e3d1a325177750eda8f6acf001642884d12264c`;
documentation fixes: `706ed2f750fef37c61594ac19de8a816382ddfac`.
Integration PR: <https://github.com/Data-Wise/probmed/pull/54>.
Review-fix PR: <https://github.com/Data-Wise/probmed/pull/55>.

## P1: default Wasserstein interval measures simulated-corner noise

`R/wasserstein-pmed.R:82-87` sets `n` to the simulated corner sample length,
then divides the influence-function variance by that number. The public
g-computation path (`R/wasserstein-pmed.R:343-354`) passes `n_mc` generated
draws to this engine and reports its SE as the data estimator's uncertainty.
It does not propagate fitted-model or empirical-covariate uncertainty.
Increasing numerical precision therefore makes the reported interval collapse
without collecting additional observations.

Reproduction on the same 500 observations, seed `20261020` for data and
`17` for each estimator call:

| n_mc | estimate | reported SE | interval width |
| ---: | ---: | ---: | ---: |
| 2,000 | 0.5009756 | 0.03780980 | 0.14821170 |
| 20,000 | 0.5720206 | 0.01181007 | 0.04629461 |
| 200,000 | 0.5597271 | 0.003726431 | 0.01460734 |

For context only, the existing row-refit bootstrap on these data, with
`n_mc=20000`, `n_boot=100`, and seed `17`, gives SE `0.07449822` and width
`0.2478665`. That small bootstrap run is not a validated oracle or a coverage
study. The blocking evidence is the public-path denominator and missing data
uncertainty, not a requirement that the two SEs agree exactly.

```r
devtools::load_all(quiet = TRUE)
set.seed(20261020)
n <- 500
C <- rnorm(n)
A <- rbinom(n, 1, plogis(0.3 * C))
M <- 0.7 * A + 0.5 * C + rnorm(n)
Y <- 0.6 * A + 0.8 * M + 0.3 * C + rnorm(n)
d <- data.frame(A, M, Y, C)
for (k in c(2000L, 20000L, 200000L)) {
  r <- wasserstein_pmed(d, n_mc = k, method = "eif", seed = 17L)
  print(c(n_mc = k, estimate = r@pmedW, se = r@pmedW_se,
          width = diff(r@pmedW_ci)))
}
wasserstein_pmed(d, n_mc = 20000L, method = "bootstrap",
                 n_boot = 100L, seed = 17L)@pmedW_se
```

The existing root-n test (`tests/testthat/test-wasserstein-pmed.R:39-52`)
validates independent, directly observed corner samples. It does not test the
public fitted-model path, where all pseudo-observations share estimated
parameters. Replacing `n_mc` by `nrow(data)` alone would not repair that
missing influence-function propagation or shared-corner dependence.

Confidence: **verified**, project scope, for the current public `method="eif"`
g-computation path away from its bootstrap fallback. This finding does not
invalidate the engine's separate direct-corner Monte Carlo calculation.
It would be invalidated by an implementation that propagates observed-data
uncertainty, with a public-path regression test and independent calibration.

## Safe next decision

Do not merge PR #54 or change its required protection yet. Choose either a
proper observed-data inference repair with regression tests and validation,
or an explicit fail-closed/experimental policy for the unvalidated EIF path.
Do not silently change the scientific default as part of a CI migration.

Review coverage: changed R implementation and CI paths were inspected, with
particular attention to single-mediator support, bootstrap, MBCO union-null
branches, numerical inversion, and new Wasserstein/RG APIs. This is not a
fresh proof of every asymptotic claim or a row-by-row audit of imported
simulation CSVs. Existing I6/I7, release, and submission holds remain.

## Validation status

Documentation fixes were regenerated with roxygen2 7.3.3; shellcheck,
actionlint, whitespace and scoped secret checks passed; six publication
scenarios passed. Full local CI was deliberately interrupted after discovery
of this blocker, so it is **not passing on either current PR head**. Its
failure statuses reflect cancellation, not an observed test assertion failure.
Main protection remains strict with `ubuntu-latest (release)` required.
