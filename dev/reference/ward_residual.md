# Gauge-Calibrated Proportion Mediated

Estimate the interventional proportion mediated together with the gauge
residual `W = R / OE` that flags non-decomposability
(treatment-by-mediator interaction). Uses a cross-fitted one-step
estimator built on the triply-robust efficient influence functions of
the four corner means `theta(a, a') = E[Y(a, M(a'))]`.

## Usage

``` r
ward_residual(
  object,
  covars = "C",
  K = 5L,
  ci_level = 0.95,
  seed = 1L,
  fieller = TRUE,
  reps = 1L,
  se_method = c("analytic", "bootstrap"),
  B = 200L,
  a0 = 0,
  a1 = 1,
  weak_id_ratio_threshold = 3,
  oe_snr_threshold = 2,
  ...
)
```

## Arguments

- object:

  A `data.frame` with columns `A` (binary treatment), `M` (mediator),
  `Y` (outcome), and the covariates named in `covars`.

- covars:

  Character vector of covariate column names. Default `"C"`.

- K:

  Integer: number of cross-fitting folds (default 5).

- ci_level:

  Numeric: confidence level (default 0.95).

- seed:

  Integer: RNG seed for fold assignment.

- fieller:

  Logical: also compute the Fieller confidence set for the ratio
  `P_med = IIE/OE`, which is unbounded when the total effect `OE` is not
  significant (default `TRUE`).

- reps:

  Integer: number of repeated cross-fitting fold draws averaged for the
  point estimate (default `1`). `reps > 1` averages the corner influence
  matrix over independent fold assignments, removing the fold-split
  component of the variance; the analytic se then adds the residual fold
  Monte-Carlo variance of the averaged point. With the corrected corner
  weights (see `se_method`) that component is small – n = 800,
  continuous: empSD 0.063 at `reps = 1` vs 0.061 at `reps = 10`, oracle
  0.056 – so the default is usually enough.

- se_method:

  Character: `"analytic"` (default) or `"bootstrap"`. Under `"analytic"`
  the se is the influence-function se of the cross-fit one-step
  estimator and the CI is symmetric Wald. Measured against exact truth
  with the corrected corner weights (`inst/sim/se_shipped_check.R`, 250
  datasets per cell, n = 800): Wald coverage 0.936 (continuous `Y`,
  strong identification), 0.972 (binary `Y`), 0.912 (intermediate `OE`),
  with the se's dispersion across datasets matching an oracle that knows
  the nuisances (CV 0.23 vs 0.18). Under `"bootstrap"` the CIs for `W`
  and `P_med` are **percentile** intervals of a nonparametric bootstrap
  that resamples rows and refits the whole cross-fit (a fresh partition
  each resample); `W_se`/`p_med` se are then bootstrap dispersion
  summaries. Re-measured with the corrected weights
  (`inst/sim/boot_gates_check.R`, B = 200, 200 datasets per cell):
  percentile coverage 0.905 (continuous) / 0.930 (binary) / 0.920
  (intermediate) at the same width as the Wald interval, and a bootstrap
  SD that is unstable for binary `Y` (CV 2.2 across datasets vs 0.33
  analytic) – an option for users who want a refit bootstrap, not a
  remedy for anything, at `B` (x `reps`) refits. The
  bootstrap-consistency result of Lin and Han (2026) for cross-fit DML
  functionals holds the nuisances fixed and does **not** cover this
  refit-per-resample scheme. `reps > 1` and `se_method = "bootstrap"`
  compose. Near the null (`OE` not bounded away from 0) the ratio is
  non-regular for every arm – even an oracle's se explodes – use the
  Fieller set.

  **History.** Every coverage figure this package reported before
  2026-08-22 – the 2,000-rep coverage grid
  (`inst/sim/results/gauge_boot_coverage_nsim2000.csv`: Wald ~0.85-0.90,
  percentile ~1.00), the 48,000-rep gate-validation grid, and the
  variance decomposition (`inst/sim/results/phi_decomp_summary.csv`;
  `docs/specs/FINDINGS-2026-08-22-phi-decomposition.md`) – was produced
  with a corner-EIF weight bug:
  [`ifelse()`](https://rdrr.io/r/base/ifelse.html) on a scalar test
  returned the **first fold row's** propensity (and mediator-density
  proxy) for every row, so the inverse-probability weights were one
  arbitrary constant per fold (NEWS, 0.3.0.9000). The "fold-split noise"
  those grids diagnosed was this bug. They are superseded. The bootstrap
  arm and both gates were re-measured with the fix the same day
  (`docs/specs/FINDINGS-2026-08-22-postfix-bootstrap-gates.md`,
  `inst/sim/results/boot_gates_postfix.csv`): see above, and the
  `weak_id` / `oe_regular` property docs of
  [GaugePmedResult](https://data-wise.github.io/probmed/dev/reference/GaugePmedResult.md).

- B:

  Integer: number of bootstrap resamples when `se_method = "bootstrap"`
  (default `200`). Cost is `B` (x `reps`) refits.

- a0, a1:

  Reference and comparison exposure levels (defaults `0`/`1`). `A` may
  use any two-level coding (factor, `{1,2}`, `{-1,1}`); it is recoded to
  the binary indicator `(A == a1)`. The gauge `W = R/OE` is invariant to
  swapping `a0`/`a1` (both `R` and `OE` flip sign). An `A` with **more
  than two** levels is an error, not a silent subset: restricting to
  `{a0,a1}` would shift the covariate-averaging population and estimate
  a different (sub-population) gauge. Filter to the two intended levels
  first. Multi-valued / continuous exposures are future work.

- weak_id_ratio_threshold:

  Numeric: the `weak_id` flag fires when the percentile `W`-CI is at
  least this many times wider than the Wald interval (default `3`).
  **Not calibrated to coverage, and cannot be** (2026-08-22, corrected
  weights): flagged and unflagged draws cover alike at every threshold
  from 1.2 to 3 – see the `weak_id` property of
  [GaugePmedResult](https://data-wise.github.io/probmed/dev/reference/GaugePmedResult.md).
  Left at 3 so the flag stays rare outside the near-null regime (0% of
  regular-cell draws, 19% near the null). Exposed so it can be tuned
  without editing the source.

- oe_snr_threshold:

  Numeric: the `oe_regular` guard flags a near-singular denominator when
  `|OE|/se(OE)` falls below this (default `2`). Validated in the same
  grid as marking **uninformatively wide** (over-covering) intervals,
  not under-coverage – see the `oe_regular` property documentation of
  [GaugePmedResult](https://data-wise.github.io/probmed/dev/reference/GaugePmedResult.md).

- ...:

  Unused.

## Value

A
[GaugePmedResult](https://data-wise.github.io/probmed/dev/reference/GaugePmedResult.md)
object.

## Examples

``` r
set.seed(1)
n <- 800; C <- rnorm(n)
A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
M <- 0.6 * A + 0.4 * C + rnorm(n)
Y <- 0.5 * A + 0.7 * M + 0.8 * A * M + 0.3 * C + rnorm(n)
ward_residual(data.frame(A, M, Y, C))
#> Gauge-calibrated proportion mediated (onestep-crossfit, n=800)
#>   P_med = 0.265  Wald [0.192, 0.338]
#>     Fieller 95% CI [0.190, 0.337]
#>   W=R/OE = 0.403  Wald [0.304, 0.502]  (p=1.22e-15)
#>   OE=1.554  IDE=0.516  IIE=0.412  R=0.626
#>   ! |W| large: additive split unreliable; interpret P_med with care.
```
