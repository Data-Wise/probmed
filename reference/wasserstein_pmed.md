# Wasserstein / Transport-Scale Proportion Mediated

Estimate `P_med^W = NIE^W / (NIE^W + NDE^W)`, the transport-scale
proportion mediated, where `NIE^W = W_2(nu_11, nu_10)` and
`NDE^W = W_2(nu_10, nu_00)` are the 2-Wasserstein distances between the
corner laws `nu(a, a')`.

Three estimation methods are available:

- `"eif"` (default): deterministic parametric g-computation with an
  observed-row influence-function SE. Includes both fitted-model and
  empirical covariate-distribution uncertainty. Requires correct working
  models, iid rows, full-rank fits, and nonzero components; d = 1 only.
  The legacy name does not assert nonparametric efficiency.

- `"bootstrap"`: row bootstrap refitting the same deterministic
  estimator. Used automatically near either zero-component boundary.
  Boundary coverage is not guaranteed by this fallback.

- `"dr"`: cross-fit DR corner-CDF estimator (multiply robust: consistent
  when either the outcome or propensity model is wrong, mediator model
  correct) with bootstrap SE. Slowest; d = 1 only.

## Usage

``` r
wasserstein_pmed(
  object,
  covars = "C",
  n_mc = 20000L,
  binY = FALSE,
  method = c("eif", "bootstrap", "dr"),
  n_boot = 200L,
  K = 5L,
  G = 300L,
  ci_level = 0.95,
  seed = 1L,
  ...
)
```

## Arguments

- object:

  A `data.frame` with columns `A` (binary 0/1), `M` (mediator), `Y`
  (outcome), and the covariates named in `covars`.

- covars:

  Character vector of covariate column names. Default `"C"`.

- n_mc:

  Legacy compatibility argument (default 20000), now ignored.
  G-computation uses deterministic quadrature instead of simulated
  corner laws.

- binY:

  Logical: if `TRUE`, treats `Y` as binary (uses logistic g-comp).

- method:

  Character: `"eif"` (default), `"bootstrap"`, or `"dr"`.

- n_boot:

  Integer: bootstrap replications (default 200). Used by `"bootstrap"`
  and `"dr"` methods, and automatically when near the boundary.

- K:

  Integer: cross-fitting folds for `method = "dr"` (default 5).

- G:

  Integer: y-grid points for `method = "dr"` (default 300).

- ci_level:

  Numeric: confidence level (default 0.95).

- seed:

  Integer: RNG seed (default 1).

- ...:

  Unused.

## Value

A
[WassersteinPmedResult](https://data-wise.github.io/probmed/reference/WassersteinPmedResult.md)
object.

## Details

The g-computation methods assume a homoscedastic Gaussian mediator model
`M ~ A + covars` and an outcome model `Y ~ A * M + covars`
(homoscedastic Gaussian, or Bernoulli-logit when `binY = TRUE`).
Covariates are standardized over the observed rows. Residual variances
use mean squared residuals, not degrees-of-freedom corrections. Gaussian
outcomes use 128-point quadrature over mixture quantiles; binary
outcomes use 64-point normal quadrature. These are numerical
approximations, not certified error bounds. Wald inference is
model-based and requires regular, nonzero distance components; it is not
a method for clustered observations or unsupported outcome families.

## Examples

``` r
set.seed(1)
n <- 800; C <- rnorm(n)
A <- rbinom(n, 1, plogis(0.3 * C))
M <- 0.7 * A + 0.5 * C + rnorm(n)
Y <- 0.6 * A + 0.8 * M + 0.3 * C + rnorm(n)
wasserstein_pmed(data.frame(A, M, Y, C), covars = "C")
#> Wasserstein P_med^W (eif, n=800)
#>   P_med^W = 0.440  [0.345, 0.535]  (se=0.048)
#>   NIE^W=0.525  NDE^W=0.668  TE^W=1.193  G^W(synergy)=0.000
```
