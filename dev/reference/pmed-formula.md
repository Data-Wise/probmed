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
  "plugin", "mbco"

- n_boot:

  Number of bootstrap samples (default: 1000)

- ci_level:

  Confidence level (default: 0.95)

- seed:

  Random seed for reproducibility

- ...:

  Additional arguments

## Value

PmedResult object

## Details

`method = "mbco"` returns a Model-Based Constrained Optimization
interval (Tofighi & Kelley, 2020): a likelihood-ratio interval for P_med
and for the indirect effect `a*b`, obtained by inverting the
constrained-likelihood test rather than by resampling. It is
deterministic (no `n_boot`, no `seed`) and supports a Gaussian outcome
and mediator, with covariates, and any treatment contrast
`x_ref != x_value`. For binary or other non-Gaussian models, use the
bootstrap methods.

For `method = "mbco"`, the `converged` flag reflects the **P_med**
interval only. The indirect-effect interval is reported separately and
may be `NA` on a degenerate design (e.g. a non-finite delta-method scale
for `a*b`) even when the P_med interval converges; check `ie_ci_lower` /
`ie_ci_upper` directly.

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
#> 95% CI: [ 0.532 ,  0.618 ]
#> 
#> Indirect Effect (Product of Coefficients):
#> Estimate: 0.147 
#> 95% CI: [ 0.0629 ,  0.248 ]
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
