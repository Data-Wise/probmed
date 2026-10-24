# Compute P_med from Formula

Compute P_med from Formula

## Arguments

- object:

  Formula for outcome model (Y ~ X + M + C)

- formula_m:

  Formula for mediator model (M ~ X + C)

- data:

  Data frame

- treatment:

  Character: treatment variable name

- mediator:

  Character: mediator variable name

- family_y:

  Family for outcome model (default: gaussian())

- family_m:

  Family for mediator model (default: gaussian())

- x_ref:

  Reference treatment value (default: 0)

- x_value:

  Treatment value (default: 1)

- method:

  Inference method: "parametric_bootstrap", "nonparametric_bootstrap",
  "plugin", "mbco", or the opt-in "model_bootstrap" (see Details).

- n_boot:

  Number of bootstrap samples (default: 1000)

- ci_level:

  Confidence level (default: 0.95)

- seed:

  Random seed for reproducibility. For single-mediator plug-in and
  bootstrap methods, an explicit seed resets the caller's RNG and the
  computation then advances it. MBCO does not use the seed.

- n_sim:

  Number of cross-world pairs drawn per evaluation of the P_med core
  (default 10000). Used by the plug-in and bootstrap methods per
  evaluation. `"model_bootstrap"` uses `100 * n_sim` pairs for its point
  estimate. Ignored by `"mbco"`, which is closed form. Distinct from
  `n_boot`.

- ...:

  Additional arguments

## Value

PmedResult object

## Details

The Gaussian MBCO exact null is the union of `a = 0` and `b = 0`; both
reduced-model likelihoods are fitted and their maximum is used.
Single-mediator MBCO uses ML residual variances (SSE/n) for both its
point estimate and likelihood. This can change the point estimate
relative to the earlier residual-df scale convention. Numerical statuses
and scan resolution are in `result@diagnostics`. Verified roots use
tolerance at most 1e-9. The bounded-domain scan has spacing at most half
the outward-search step and can miss narrower components; it does not
prove global connectedness. Detected disconnected sets or failed
profiles return missing interval bounds. Domain-truncated bounds remain
diagnostic bounds and are not verified LR endpoints.

`method = "model_bootstrap"` (opt-in; not the default) simulates
complete synthetic datasets from the fitted sequential models with the
observed treatment and covariate rows held fixed: the mediator from its
model, then the outcome using the simulated mediator. Both equations are
refitted on the same rows, so coefficients and residual scales are
re-estimated jointly. The estimate is the fitted-model plug-in; the
percentile interval is reported only when at least 99% of draws are
valid (finite, converged, not separated), otherwise both bounds are `NA`
with `result@diagnostics$interval_status = "bootstrap_incomplete"`. It
models uncertainty under the specified model; it does not correct
misspecification. Calibration is limited to the documented simulation
cells (see NEWS); it does not establish general nominal coverage.

`method = "mbco"` returns a Model-Based Constrained Optimization
interval (Tofighi & Kelley, 2020): a likelihood-ratio interval for P_med
and for the indirect effect `a*b`, obtained by inverting the
constrained-likelihood test rather than by resampling. It is
deterministic (no `n_boot`, no `seed`) and supports a Gaussian outcome
and mediator, with covariates, and any treatment contrast
`x_ref != x_value`. For binary or other non-Gaussian models, use the
bootstrap methods.

For single-mediator `method = "mbco"`, `converged` requires two verified
**P_med** LR endpoints and no detected re-entry on the scan. It does not
establish coverage calibration. The indirect-effect interval is reported
separately and may be `NA` on a degenerate design (e.g. a non-finite
delta-method scale for `a*b`) even when the P_med interval converges;
check `ie_ci_lower` / `ie_ci_upper` directly.

**Composite null.** No mediation is the union of a = 0 and b = 0. At the
double null a = b = 0, the profile LR at P_med = 0.5 is the minimum of
the two branch statistics, so the chi-square(1) interval is
conservative: calibration coverage was 0.997 to 0.998 at the 95% level,
and about 0.994 near the origin. `@diagnostics$pmed$branch_lr` holds the
LR statistics for a = 0 and b = 0. `@diagnostics$pmed$near_double_null`
is `TRUE` when both are below the cutoff (neither path individually
significant), and [`print()`](https://rdrr.io/r/base/print.html) then
notes that the interval is conservative. Minimum observed coverage on
the tested lattice was 0.9415, not a uniform coverage guarantee. It
stays conservative until one path is about 4 standard errors from zero,
so a single null with a weak other path is conservative too, and the
flag can be `FALSE` there.

## Examples

``` r
# Toy example: Simple mediation model
# Generate data where X affects Y through M
set.seed(123)
n <- 100
data <- data.frame(
  X = rnorm(n),
  C = rnorm(n)
)
data$M <- 0.5 * data$X + 0.3 * data$C + rnorm(n, sd = 0.5)
data$Y <- 0.4 * data$M + 0.2 * data$X + 0.2 * data$C + rnorm(n, sd = 0.5)

# Compute P_med using plugin estimator (fast, no CI)
result_plugin <- pmed(
  Y ~ X + M + C,
  formula_m = M ~ X + C,
  data = data,
  treatment = "X",
  mediator = "M",
  method = "plugin"
)
print(result_plugin)
#> 
#> P_med: Probability of Mediated Shift
#> ====================================
#> 
#> Estimate: 0.574 
#> 
#> Indirect Effect (Product of Coefficients):
#> Estimate: 0.148 
#> 
#> Inference: plugin 
#> 
#> Treatment contrast: X = 1 vs. X* = 0 
#> 
#> Interpretation:
#>   P(Y(1, M(1)) > Y(1, M(0))) = 0.574
#>   P that the mediator shift (M(0) -> M(1)) leaves a random individual better off, holding X = 1.
#> 

# \donttest{
# With parametric bootstrap for confidence intervals
result_boot <- pmed(
  Y ~ X + M + C,
  formula_m = M ~ X + C,
  data = data,
  treatment = "X",
  mediator = "M",
  method = "parametric_bootstrap",
  n_boot = 200,  # Use more (e.g., 1000+) in practice
  seed = 456
)
print(result_boot)
#> 
#> P_med: Probability of Mediated Shift
#> ====================================
#> 
#> Estimate: 0.574 
#> 95% CI: [ 0.533 ,  0.618 ]
#> 
#> Indirect Effect (Product of Coefficients):
#> Estimate: 0.147 
#> 95% CI: [ 0.065 ,  0.248 ]
#> 
#> Inference: parametric_bootstrap 
#> Bootstrap samples: 200 
#> 
#> Treatment contrast: X = 1 vs. X* = 0 
#> 
#> Interpretation:
#>   P(Y(1, M(1)) > Y(1, M(0))) = 0.574
#>   P that the mediator shift (M(0) -> M(1)) leaves a random individual better off, holding X = 1.
#> 
# }
```
