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
  "plugin"

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
#> Estimate: 0.391 
#> 
#> Inference: plugin 
#> 
#> Treatment contrast: X = 1 vs. X* = 0 
#> 
#> Interpretation:
#>   P(Y_{X*, M_X} > Y_{X, M_X}) = 0.391 
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
#> Estimate: 0.277 
#> 95% CI: [ 0.231 ,  0.329 ]
#> 
#> Inference: parametric_bootstrap 
#> Bootstrap samples: 200 
#> 
#> Treatment contrast: X = 1 vs. X* = 0 
#> 
#> Interpretation:
#>   P(Y_{X*, M_X} > Y_{X, M_X}) = 0.277 
#> 
# }
```
