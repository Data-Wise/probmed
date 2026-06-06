# Working with lavaan Models

The `probmed` package integrates seamlessly with `lavaan`, the popular R
package for structural equation modeling (SEM). This vignette
demonstrates how to extract mediation structures from `lavaan` models
and compute the probabilistic index of mediation ($`P_{med}`$).

## Prerequisites

Ensure you have both `probmed` and `lavaan` installed:

``` r

install.packages("lavaan")
# install.packages("probmed") # Not on CRAN yet
```

Load the packages:

``` r

library(probmed)
library(lavaan)
```

## Simple Mediation Model

Let’s start with a simple mediation model where $`X`$ affects $`M`$, and
both $`X`$ and $`M`$ affect $`Y`$.

### 1. Simulate Data

First, we generate some data consistent with a mediation process:

``` r

set.seed(123)
n <- 200
X <- rnorm(n)
M <- 0.5 * X + rnorm(n)
Y <- 0.4 * M + 0.2 * X + rnorm(n)
data <- data.frame(X = X, M = M, Y = Y)
```

### 2. Fit Model with lavaan

We specify the model using `lavaan` syntax. Note that we don’t strictly
need to label the parameters, but it’s good practice.

``` r

model <- '
  # Mediator model
  M ~ a*X
  
  # Outcome model
  Y ~ b*M + cp*X
  
  # Variances (optional, lavaan adds them by default)
  M ~~ M
  Y ~~ Y
'

fit <- sem(model, data = data)
summary(fit)
```

### 3. Extract Mediation Structure

Use
[`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md)
to convert the `lavaan` object into a format `probmed` understands. You
need to specify the treatment and mediator variables. The outcome
variable is usually auto-detected, but you can specify it explicitly if
needed.

``` r

extract <- extract_mediation(fit, treatment = "X", mediator = "M")
print(extract)
```

### 4. Compute $`P_{med}`$

Now you can compute $`P_{med}`$ using any of the supported methods
(plugin, parametric bootstrap, or non-parametric bootstrap).

``` r

# Plugin estimator (fast, point estimate only)
pmed(extract, method = "plugin")

# Parametric bootstrap (with confidence intervals)
# Note: Increase n_boot for production use
pmed(extract, method = "parametric_bootstrap", n_boot = 200, seed = 456)
```

## Models with Covariates

`probmed` automatically handles covariates included in the `lavaan`
model.

``` r

# Add covariates
C1 <- rnorm(n)
C2 <- rnorm(n)
data_cov <- data.frame(X = X, M = M, Y = Y, C1 = C1, C2 = C2)

# Model with covariates
model_cov <- '
  M ~ X + C1 + C2
  Y ~ M + X + C1 + C2
'

fit_cov <- sem(model_cov, data = data_cov)

# Extract
extract_cov <- extract_mediation(fit_cov, treatment = "X", mediator = "M")

# Compute P_med
pmed(extract_cov, method = "plugin")
```

## Standardized Estimates

If you prefer to use standardized estimates, you can set
`standardized = TRUE` in
[`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md).
Note that this affects the interpretation of the path coefficients.

``` r

extract_std <- extract_mediation(
  fit, 
  treatment = "X", 
  mediator = "M", 
  standardized = TRUE
)

print(extract_std)
```

## Troubleshooting

### “Outcome variable not found”

If `probmed` cannot auto-detect the outcome variable (e.g., if there are
multiple endogenous variables besides the mediator), you must specify it
explicitly:

``` r

extract_mediation(fit, treatment = "X", mediator = "M", outcome = "Y")
```

### “Path not found”

Ensure that the model actually contains the paths $`M \sim X`$,
$`Y \sim M`$, and $`Y \sim X`$. If any of these are missing (e.g., if
you fixed $`c'`$ to 0),
[`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md)
will fail because $`P_{med}`$ requires all three paths (though $`c'`$
can be 0, it must be in the model structure or handled as 0).

Currently, `probmed` expects the direct effect $`Y \sim X`$ to be
present in the model. If you omitted it (assuming full mediation), you
can add it and fix it to 0 in `lavaan`: `Y ~ 0*X`.
