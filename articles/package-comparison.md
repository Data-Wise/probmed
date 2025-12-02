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
#> <probmed::PmedResult>
#>  @ estimate         : num 0.483
#>  @ ci_lower         : Named num 0.417
#>  .. - attr(*, "names")= chr "2.5%"
#>  @ ci_upper         : Named num 0.551
#>  .. - attr(*, "names")= chr "97.5%"
#>  @ ci_level         : num 0.95
#>  @ method           : chr "parametric_bootstrap"
#>  @ n_boot           : int 1000
#>  @ boot_estimates   : num [1:1000] 0.479 0.48 0.506 0.513 0.5 ...
#>  @ ie_estimate      : num 0.189
#>  @ ie_ci_lower      : Named num 0.0809
#>  .. - attr(*, "names")= chr "2.5%"
#>  @ ie_ci_upper      : Named num 0.326
#>  .. - attr(*, "names")= chr "97.5%"
#>  @ ie_boot_estimates: num [1:1000] 0.219 0.208 0.169 0.269 0.225 ...
#>  @ x_ref            : num 0
#>  @ x_value          : num 1
#>  @ source_extract   : <probmed::MediationExtract>
#>  .. @ estimates          : Named num [1:5] -0.1028 0.4475 0.1351 0.4238 0.0549
#>  .. .. - attr(*, "names")= chr [1:5] "(Intercept)" "X" "(Intercept)" "M" ...
#>  .. @ mediator_predictors: chr [1:2] "(Intercept)" "X"
#>  .. @ outcome_predictors : chr [1:3] "(Intercept)" "M" "X"
#>  .. @ a_path             : num 0.448
#>  .. @ b_path             : num 0.424
#>  .. @ c_prime            : num 0.0549
#>  .. @ vcov               : num [1:5, 1:5] 0.00952 -0.00103 0 0 0 ...
#>  .. @ sigma_m            : num 0.971
#>  .. @ sigma_y            : num 0.951
#>  .. @ data               :'data.frame':  100 obs. of  2 variables:
#>  .. .. $ M: num  -0.991 0.142 0.533 -0.312 -0.887 ...
#>  .. .. $ X: num  -0.5605 -0.2302 1.5587 0.0705 0.1293 ...
#>  .. .. - attr(*, "terms")=Classes 'terms', 'formula'  language M ~ X
#>  .. ..  .. ..- attr(*, "variables")= language list(M, X)
#>  .. ..  .. ..- attr(*, "factors")= int [1:2, 1] 0 1
#>  .. ..  .. .. ..- attr(*, "dimnames")=List of 2
#>  .. ..  .. .. .. ..$ : chr [1:2] "M" "X"
#>  .. ..  .. .. .. ..$ : chr "X"
#>  .. ..  .. ..- attr(*, "term.labels")= chr "X"
#>  .. ..  .. ..- attr(*, "order")= int 1
#>  .. ..  .. ..- attr(*, "intercept")= int 1
#>  .. ..  .. ..- attr(*, "response")= int 1
#>  .. ..  .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
#>  .. ..  .. ..- attr(*, "predvars")= language list(M, X)
#>  .. ..  .. ..- attr(*, "dataClasses")= Named chr [1:2] "numeric" "numeric"
#>  .. ..  .. .. ..- attr(*, "names")= chr [1:2] "M" "X"
#>  .. @ n_obs              : int 100
#>  .. @ source_package     : chr "stats"
#>  .. @ converged          : logi TRUE
#>  .. @ treatment          : chr "X"
#>  .. @ mediator           : chr "M"
#>  .. @ outcome            : chr "Y"
#>  @ converged        : logi TRUE
#>  @ call             : NULL
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
