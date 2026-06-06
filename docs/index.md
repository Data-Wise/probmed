# probmed

**probmed** provides a robust framework for computing **$`P_{med}`$**, a
scale-free probabilistic effect size for causal mediation analysis.

## What is $`P_{med}`$?

$`P_{med}`$ offers an intuitive alternative to traditional mediation
effect sizes. It represents the probability that the outcome for an
individual in the treatment group is higher than the outcome for an
individual in the control group, specifically due to the indirect path
through the mediator.

Unlike the “Proportion Mediated” or simple indirect effect coefficients,
$`P_{med}`$: \* **Is Scale-Free**: Interpretable regardless of the units
of measurement for $`X`$, $`M`$, or $`Y`$. \* **Is Bounded**: Always
falls between 0 and 1 (or 0 and 0.5 for some definitions), making it
comparable across studies. \* **Handles Non-Linearity**: Can be extended
to Generalized Linear Models (GLMs) where coefficients are not directly
comparable.

## Features

- **Flexible Modeling**: Supports both **Linear Models (LM)** and
  **Generalized Linear Models (GLM)** (e.g., Logistic, Poisson) for the
  mediator and outcome.
- **Comprehensive Effect Reporting**: Provides both $`P_{med}`$
  (probabilistic effect size) and the traditional Indirect Effect
  ($`a \times b`$) with bootstrap confidence intervals.
- **Robust Inference**: Implements multiple methods for constructing
  confidence intervals:
  - **Parametric Bootstrap**: Efficient and accurate for large samples
    with known distributions.
  - **Nonparametric Bootstrap**: Robust to distributional assumptions
    (resampling-based).
  - **Plug-in Estimator**: Fast point estimation.
- **Modern Architecture**: Built on the **S7** object-oriented system,
  ensuring type safety, stability, and easy extensibility.
- **Ecosystem Integration**: Works seamlessly with popular R packages:
  - **lavaan**: Extract from Structural Equation Models with FIML and
    robust estimators
  - **mediation**: Direct integration with
    [`mediate()`](https://rdrr.io/pkg/mediation/man/mediate.html)
    objects
  - Standard **lm/glm**: Native support for regression models
- **User-Friendly Interface**: Designed to work with standard R
  `formula` interfaces and `data.frame` inputs.

## Installation

You can install the development version of probmed from
[GitHub](https://github.com/) with:

``` r

# install.packages("devtools")
devtools::install_github("data-wise/probmed")
```

## Quick Example

Here is a basic example showing how to simulate data and compute
$`P_{med}`$ using a parametric bootstrap approach.

### 1. Simulate Data

First, we generate a dataset with a continuous treatment ($`X`$),
mediator ($`M`$), and outcome ($`Y`$), along with a covariate ($`C`$).

``` r

library(probmed)

set.seed(123)
n <- 300
data <- data.frame(X = rnorm(n), C = rnorm(n))

# Mediator model: M ~ X + C
data$M <- 0.5 * data$X + 0.3 * data$C + rnorm(n)

# Outcome model: Y ~ M + X + C
data$Y <- 0.4 * data$M + 0.2 * data$X + 0.2 * data$C + rnorm(n)
```

### 2. Compute $`P_{med}`$

We use the
[`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md)
function to estimate the effect size. We specify the models for the
outcome and mediator, and select the bootstrap method.

``` r

result <- pmed(
  Y ~ X + M + C,             # Outcome model formula
  formula_m = M ~ X + C,     # Mediator model formula
  data = data,
  treatment = "X",
  mediator = "M",
  method = "parametric_bootstrap",
  n_boot = 1000,
  seed = 42
)
```

### 3. View Results

The results object provides a clear summary of both $`P_{med}`$ and the
Indirect Effect with confidence intervals.

``` r

print(result)
#>
#> P_med: Probability of Mediated Shift
#> ====================================
#>
#> Estimate: 0.563
#> 95% CI: [0.520, 0.605]
#>
#> Indirect Effect (Product of Coefficients):
#> Estimate: 0.198
#> 95% CI: [0.134, 0.268]
#>
#> Inference: parametric_bootstrap
#> Bootstrap samples: 1000
#>
#> Treatment contrast: X = 1 vs. X* = 0
#>
#> Interpretation:
#>   P(Y_{X*, M_X} > Y_{X, M_X}) = 0.563

summary(result)
```

**Interpretation**: - **$`P_{med}`$ = 0.563**: There is a 56.3%
probability that a treated individual will have a higher outcome than a
control individual through the indirect path. Since 0.50 represents “no
effect” (random chance), this indicates a positive mediation effect. -
**Indirect Effect = 0.198**: The traditional product-of-coefficients
indirect effect ($`a \times b`$) quantifies the average change in
outcome due to the mediated path.

## Workflow Options

`probmed` offers three workflow approaches to suit your needs:

### 1. Formula Interface (Recommended for Most Users)

Directly specify your models using R formulas:

``` r

result <- pmed(
  Y ~ X + M + C,             # Outcome model
  formula_m = M ~ X + C,     # Mediator model
  data = data,
  treatment = "X",
  mediator = "M",
  method = "parametric_bootstrap"
)
```

### 2. Integration with lavaan (For SEM Users)

Extract from fitted Structural Equation Models with full support for
FIML, robust estimators, and standardized estimates:

``` r

library(lavaan)

# Define SEM model
model <- '
  M ~ a*X + C
  Y ~ b*M + cp*X + C
'

# Fit with FIML for missing data
fit <- sem(model, data = data, missing = "fiml")

# Extract and compute P_med
extract <- extract_mediation(fit, treatment = "X", mediator = "M")
result <- pmed(extract, method = "parametric_bootstrap", n_boot = 5000)
print(result)
```

### 3. Integration with mediation Package

Seamlessly work with objects from the `mediation` package:

``` r

library(mediation)

# Fit models
model_m <- lm(M ~ X + C, data = data)
model_y <- lm(Y ~ M + X + C, data = data)

# Run mediate
med_out <- mediate(model_m, model_y, treat = "X", mediator = "M", boot = TRUE)

# Compute P_med from mediate object
result <- pmed(extract_mediation(med_out))
print(result)
```

## Advanced Features

### Generalized Linear Models (GLMs)

`probmed` supports GLMs for non-normal outcomes. For example, with a
binary outcome:

``` r

# Binary outcome example
data_bin <- data
data_bin$Y_binary <- rbinom(n, 1, plogis(0.4 * data$M + 0.2 * data$X))

result_bin <- pmed(
  Y_binary ~ X + M + C,
  formula_m = M ~ X + C,
  data = data_bin,
  family_y = binomial(),    # Logistic regression for binary outcome
  treatment = "X",
  mediator = "M",
  method = "parametric_bootstrap",
  n_boot = 1000
)
print(result_bin)
```

### Bootstrap Methods Comparison

**When to use each method:**

- **Plugin** (`method = "plugin"`): Fast point estimates only, no
  confidence intervals. Use for quick exploration.
- **Parametric Bootstrap** (`method = "parametric_bootstrap"`): Default
  choice. Fast and accurate when parameter estimates are approximately
  normally distributed.
- **Nonparametric Bootstrap** (`method = "nonparametric_bootstrap"`):
  More robust to distributional assumptions but computationally
  intensive. Use when parametric assumptions are questionable.

``` r

# Compare methods
result_plugin <- pmed(Y ~ X + M, formula_m = M ~ X, data = data,
                      treatment = "X", mediator = "M", method = "plugin")

result_param <- pmed(Y ~ X + M, formula_m = M ~ X, data = data,
                     treatment = "X", mediator = "M",
                     method = "parametric_bootstrap", n_boot = 1000)

result_nonparam <- pmed(Y ~ X + M, formula_m = M ~ X, data = data,
                        treatment = "X", mediator = "M",
                        method = "nonparametric_bootstrap", n_boot = 1000)
```

## References

- Tofighi, D. (In Press). *Probabilistic Effect Sizes for Mediation
  Analysis*.
- MacKinnon, D. P. (2008). *Introduction to Statistical Mediation
  Analysis*. Lawrence Erlbaum Associates.

## Contributing

We welcome contributions! Please see the [Issue
Tracker](https://github.com/data-wise/probmed/issues) for potential
features or bug reports.
