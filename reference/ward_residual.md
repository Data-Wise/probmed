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
  Monte-Carlo variance of the averaged point.

- se_method:

  Character: `"analytic"` (default, influence-function se, symmetric
  Wald CI) or `"bootstrap"`. The analytic se for the ratios `W` and
  `P_med` is right-skewed with a median below the empirical SD, so the
  symmetric Wald CI is mildly **anti-conservative** (sub-nominal
  coverage \\\approx 0.85\\-\\0.90\\). Under `"bootstrap"` the CIs for
  `W` and `P_med` are the tail-aware **percentile** intervals of the
  nonparametric bootstrap (resample rows, refit) – the appropriate
  construction for a skewed ratio (widening a symmetric se fails to
  cover it); `W_se`/`p_med` se are still reported as bootstrap
  dispersion summaries. `reps > 1` and `se_method = "bootstrap"`
  compose. Bootstrap validity follows from the estimator being a
  Neyman-orthogonal cross-fit (DML-type) functional (Lin et al. 2026);
  it requires the ratio to be regular, i.e. `OE` bounded away from 0 –
  see the Fieller diagnostic for the near-null case.

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

- ...:

  Unused.

## Value

A
[GaugePmedResult](https://data-wise.github.io/probmed/reference/GaugePmedResult.md)
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
#>   P_med = 0.271  Wald [0.210, 0.332]
#>     Fieller 95% CI [0.206, 0.330]
#>   W=R/OE = 0.375  Wald [0.190, 0.561]  (p=7.44e-05)
#>   OE=1.566  IDE=0.554  IIE=0.424  R=0.588
#>   ! |W| large: additive split unreliable; interpret P_med with care.
```
