# Comparing probmed Workflows

The `probmed` package is designed to be flexible, allowing you to
extract mediation structures from various R objects. This vignette
compares the different workflows available:

1.  **Base R (`lm`/`glm`)**: Best for simple models and maximum control.
2.  **`lavaan`**: Best for Structural Equation Models (SEM), latent
    variables, and complex path structures.
3.  **`mediation`**: Best if you are already using the `mediation`
    package for other analyses (e.g., sensitivity analysis).

## Feature Comparison

| Feature                | `lm` / `glm` |           `lavaan`           |      `mediation`       |
|:-----------------------|:------------:|:----------------------------:|:----------------------:|
| **Model Types**        | Linear, GLM  |   SEM, CFA, Path Analysis    | Linear, GLM, GAM, etc. |
| **Latent Variables**   |      No      |             Yes              |           No           |
| **Multiple Mediators** | Manual setup |        Native support        |     Native support     |
| **Missing Data**       | `na.action`  |   FIML (`missing="fiml"`)    |      `na.action`       |
| **Standardized Est.**  |    Manual    | Native (`standardized=TRUE`) |         Manual         |
| **Ease of Use**        |     High     |            Medium            |          High          |

## Workflow 1: Base R (`lm`/`glm`)

This is the most direct method. You fit two separate models: one for the
mediator and one for the outcome.

``` r
library(probmed)

# Simulate data
set.seed(123)
n <- 100
data <- data.frame(X = rnorm(n))
data$M <- 0.5 * data$X + rnorm(n)
data$Y <- 0.4 * data$M + 0.2 * data$X + rnorm(n)

# Fit models
model_m <- lm(M ~ X, data = data)
model_y <- lm(Y ~ M + X, data = data)

# Extract and compute
extract <- extract_mediation(model_m, model_y = model_y, treatment = "X", mediator = "M")
pmed(extract)
```

**Pros**: \* No extra dependencies. \* Full control over model
specification.

**Cons**: \* Requires managing two model objects. \* Manual handling of
missing data consistency.

## Workflow 2: `lavaan`

Ideal for users familiar with SEM or when dealing with latent variables
or FIML.

``` r
library(lavaan)

model <- '
  M ~ a*X
  Y ~ b*M + cp*X
'

fit <- sem(model, data = data)

extract <- extract_mediation(fit, treatment = "X", mediator = "M")
pmed(extract)
```

**Pros**: \* Single model specification. \* Handles FIML for missing
data. \* Supports latent variables (though `probmed` currently requires
observed data for bootstrapping).

**Cons**: \* Requires `lavaan` package. \* Syntax is different from
standard R formulas.

## Workflow 3: `mediation`

If you use the `mediation` package for its specific features (like
sensitivity analysis), you can easily pipe the results to `probmed`.

``` r
library(mediation)

med_out <- mediate(model_m, model_y, treat = "X", mediator = "M")

extract <- extract_mediation(med_out)
pmed(extract)
```

**Pros**: \* Seamless integration with `mediation` workflow. \*
Automatic variable detection.

**Cons**: \* Requires `mediation` package. \* Limited to models
supported by
[`mediate()`](https://rdrr.io/pkg/mediation/man/mediate.html).

## Conclusion

Choose the workflow that best fits your existing analysis pipeline.
`probmed` ensures that regardless of the input method, the calculation
of $P_{med}$ remains consistent and robust.
