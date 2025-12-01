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
if (FALSE) { # \dontrun{
# Simulate data
set.seed(123)
n <- 300
data <- data.frame(
  X = rnorm(n),
  C = rnorm(n)
)
data$M <- 0.5 * data$X + 0.3 * data$C + rnorm(n)
data$Y <- 0.4 * data$M + 0.2 * data$X + 0.2 * data$C + rnorm(n)

# Compute P_med
result <- pmed(
  Y ~ X + M + C,
  formula_m = M ~ X + C,
  data = data,
  treatment = "X",
  mediator = "M"
)
print(result)
} # }
```
