# Sobol / Variance-Scale Proportion Mediated

Estimate the **Sobol proportion mediated**
`P_med^{sigma2} = V_med / V_T` – the share of the total interventional
outcome variance carried by the mediated pathway, i.e. the first-order
Sobol (functional-ANOVA) sensitivity index of the outcome with respect
to the mediator. The estimator is cross-fitted and built on the same
triply-robust corner pseudo-outcomes as
[`ward_residual()`](https://data-wise.github.io/probmed/reference/ward_residual.md)
and
[`incr_pmed()`](https://data-wise.github.io/probmed/reference/incr_pmed.md);
the standard error uses the delta method through the quadratic variance
components and the ratio identity. See
[SobolPmedResult](https://data-wise.github.io/probmed/reference/SobolPmedResult.md)
for the variance decomposition and the boundary remedy.

## Usage

``` r
sobol_pmed(
  object,
  pd = 0.5,
  pm = 0.5,
  covars = "C",
  K = 5L,
  seed = 1L,
  ci_level = 0.95,
  warn_boundary = TRUE,
  boundary_test = c("split", "plugin", "none"),
  procedure = c("B", "A"),
  reps = 1L,
  se_method = c("analytic", "bootstrap"),
  B = 200L,
  ...
)
```

## Arguments

- object:

  A `data.frame` with columns `A` (binary treatment), `M` (mediator),
  `Y` (continuous outcome), and the covariates named in `covars`.

- pd, pm:

  Numeric: tilt probabilities for the direct and mediated pathways
  (Bernoulli design points); defaults `0.5`.

- covars:

  Character vector of covariate column names. Default `"C"`.

- K:

  Integer: number of cross-fitting folds (default `5`).

- seed:

  Integer: RNG seed for fold assignment and the boundary split (default
  `1`).

- ci_level:

  Numeric: confidence level (default `0.95`).

- warn_boundary:

  Logical: warn when the boundary test does not reject and the reported
  interval falls back to the Procedure-A upper bound (default `TRUE`).

- boundary_test:

  Character: boundary test for `H0: V_med = 0`. `"split"` (default) is
  the Williamson sample-split test on `Delta_m`; `"plugin"` is the
  (over-rejecting) plug-in Wald test, diagnostic only.

- procedure:

  Character: interval procedure for `ci`. `"B"` (default) is the
  sample-split CI (image of the `Delta_m` CI); `"A"` is the legacy gated
  interval (non-uniform near the boundary). See
  [SobolPmedResult](https://data-wise.github.io/probmed/reference/SobolPmedResult.md).

- reps:

  Integer: repeated cross-fitting draws averaged for the point estimate
  (default `1`). Use `reps > 1` near the boundary for a reproducible,
  lower-variance estimate (it removes the ~80% fold-split variance).

- se_method:

  Character: `"analytic"` (default) or `"bootstrap"`. The bootstrap
  gives a valid (conservative) se near the non-regular boundary, where
  the analytic se is anti-conservative; off the boundary the analytic se
  is calibrated.

- B:

  Integer: bootstrap resamples when `se_method = "bootstrap"` (default
  `200`).

- ...:

  Unused.

## Value

A
[SobolPmedResult](https://data-wise.github.io/probmed/reference/SobolPmedResult.md)
object.

## Details

At the null `V_med = 0` (equivalently `Delta_m = 0`) the
squared-variance influence function degenerates and the symmetric Wald
interval is non-regular. A Williamson et al. (2021) sample-split test on
the regular contrast `Delta_m` decides the boundary. The default
reported interval is **Procedure B** (`procedure = "B"`): the image of
the regular `Delta_m` Wald CI under `delta -> (c_m/V_T) delta^2`
(continuous mapping, no gate). Its validity inherits that of the input
`Delta_m` CI – so near the boundary use `reps > 1` and
`se_method = "bootstrap"` for a valid interval (the analytic `Delta_m`
se is anti-conservative there; see the boundary coverage section).
**Procedure A** (`procedure = "A"`) is the legacy gated rule (Wald on
rejection, one-sided `[0, Pmed_upper]` otherwise), pointwise- but
**not** uniformly-valid across the near-null transition (Leeb-Potscher
pre-test under-coverage). The boundary machinery recurses through a
plain internal fitter (not this generic), so the point estimate, fold
draw, and split are deterministic in `seed`. See the boundary coverage
section in
[SobolPmedResult](https://data-wise.github.io/probmed/reference/SobolPmedResult.md).

## References

Williamson, B. D., Gilbert, P. B., Carone, M., & Simon, N. (2021).
Nonparametric variable importance assessment using machine learning
techniques. *Biometrics*, 77(1), 9–22.

## Examples

``` r
set.seed(1)
n <- 1500; C <- rnorm(n)
A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
M <- 0.6 * A + 0.4 * C + rnorm(n)
# strong mediation (off-boundary): Wald interval reported
Y <- 0.5 * A + 0.7 * M + 0.3 * C + rnorm(n)
sobol_pmed(data.frame(A, M, Y, C))
#> Sobol / variance-scale proportion mediated P_med^sigma2 (onestep-crossfit, n=1500)
#>   P_med = 0.407  se = 0.104
#>   Wald 95% CI [0.221, 0.649]
#>   H0 V_med=0 rejected (split p=0.0114).
#>   S1_med = 0.407  ST_med = 0.416
#>   Vd=0.0503  Vm=0.0350  Vint=0.0008  VT=0.0861
```
