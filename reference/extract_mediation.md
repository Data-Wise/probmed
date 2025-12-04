# Extract Mediation Structure

Generic function to extract mediation structure (a, b, c' paths and
variance-covariance matrices) from fitted models. Supports multiple
model types including lm/glm, lavaan SEM objects, and mediation package
objects.

## Usage

``` r
extract_mediation(object, ...)
```

## Arguments

- object:

  Fitted model object (lm, glm, lavaan, or mediate)

- ...:

  Additional arguments passed to methods. Common arguments include:

  - `treatment`: Character name of treatment variable

  - `mediator`: Character name of mediator variable

  - `outcome`: Character name of outcome variable (auto-detected for
    lavaan)

  - `model_y`: Outcome model (for lm/glm methods)

  - `data`: Data frame (for lm/glm methods)

## Value

`MediationExtract` object containing:

- Path coefficients (a, b, c')

- Variance-covariance matrix

- Residual standard deviations

- Variable names and metadata

## See also

[`pmed`](https://data-wise.github.io/probmed/reference/pmed.md) for
computing P_med from extracted objects

## Examples

``` r
# Example 1: Extract from lm objects
set.seed(123)
data <- data.frame(X = rnorm(100), C = rnorm(100))
data$M <- 0.5 * data$X + 0.3 * data$C + rnorm(100)
data$Y <- 0.4 * data$M + 0.2 * data$X + rnorm(100)

model_m <- lm(M ~ X + C, data = data)
model_y <- lm(Y ~ M + X + C, data = data)
extract <- extract_mediation(model_m, model_y = model_y,
                              treatment = "X", mediator = "M", data = data)
print(extract)
#> <probmed::MediationExtract>
#>  @ estimates          : Named num [1:7] 0.1351 0.3668 0.3238 -0.0193 0.3426 ...
#>  .. - attr(*, "names")= chr [1:7] "(Intercept)" "X" "C" "(Intercept)" ...
#>  @ mediator_predictors: chr [1:3] "(Intercept)" "X" "C"
#>  @ outcome_predictors : chr [1:4] "(Intercept)" "M" "X" "C"
#>  @ a_path             : num 0.367
#>  @ b_path             : num 0.343
#>  @ c_prime            : num 0.173
#>  @ vcov               : num [1:7, 1:7] 0.009243 -0.000939 0.001007 0 0 ...
#>  @ sigma_m            : num 0.951
#>  @ sigma_y            : num 1.05
#>  @ data               :'data.frame': 100 obs. of  4 variables:
#>  .. $ X: num  -0.5605 -0.2302 1.5587 0.0705 0.1293 ...
#>  .. $ C: num  -0.71 0.257 -0.247 -0.348 -0.952 ...
#>  .. $ M: num  1.705 1.274 0.44 0.474 -0.635 ...
#>  .. $ Y: num  -0.145 -0.289 -0.451 -0.849 -0.665 ...
#>  @ n_obs              : int 100
#>  @ source_package     : chr "stats"
#>  @ converged          : logi TRUE
#>  @ treatment          : chr "X"
#>  @ mediator           : chr "M"
#>  @ outcome            : chr "Y"

# \donttest{
# Example 2: Extract from lavaan
if (requireNamespace("lavaan", quietly = TRUE)) {
  library(lavaan)
  model <- 'M ~ a*X + C
            Y ~ b*M + cp*X + C'
  fit <- sem(model, data = data)
  extract <- extract_mediation(fit, treatment = "X", mediator = "M")
  print(extract)
}
#> This is lavaan 0.6-20
#> lavaan is FREE software! Please report any bugs.
#> <probmed::MediationExtract>
#>  @ estimates          : Named num [1:3] 0.367 0.343 0.173
#>  .. - attr(*, "names")= chr [1:3] "a" "b" "cp"
#>  @ mediator_predictors: chr [1:2] "X" "C"
#>  @ outcome_predictors : chr [1:3] "M" "X" "C"
#>  @ a_path             : num 0.367
#>  @ b_path             : num 0.343
#>  @ c_prime            : num 0.173
#>  @ vcov               : num [1:3, 1:3] 1.07e-02 3.77e-18 -5.21e-18 3.65e-18 1.21e-02 ...
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:3] "a" "b" "cp"
#>  ..  ..$ : chr [1:3] "a" "b" "cp"
#>  @ sigma_m            : num 0.937
#>  @ sigma_y            : num 1.03
#>  @ data               :'data.frame': 100 obs. of  4 variables:
#>  .. $ M: num  1.705 1.274 0.44 0.474 -0.635 ...
#>  .. $ Y: num  -0.145 -0.289 -0.451 -0.849 -0.665 ...
#>  .. $ X: num  -0.5605 -0.2302 1.5587 0.0705 0.1293 ...
#>  .. $ C: num  -0.71 0.257 -0.247 -0.348 -0.952 ...
#>  @ n_obs              : int 100
#>  @ source_package     : chr "lavaan"
#>  @ converged          : logi TRUE
#>  @ treatment          : chr "X"
#>  @ mediator           : chr "M"
#>  @ outcome            : chr "Y"

# Example 3: Extract from mediate object
if (requireNamespace("mediation", quietly = TRUE)) {
  library(mediation)
  med_out <- mediation::mediate(model_m, model_y,
                                 treat = "X", mediator = "M", boot = FALSE)
  extract <- extract_mediation(med_out)
  print(extract)
}
#> Loading required package: MASS
#> Loading required package: Matrix
#> Loading required package: mvtnorm
#> Loading required package: sandwich
#> mediation: Causal Mediation Analysis
#> Version: 4.5.1
#> <probmed::MediationExtract>
#>  @ estimates          : Named num [1:7] 0.1351 0.3668 0.3238 -0.0193 0.3426 ...
#>  .. - attr(*, "names")= chr [1:7] "(Intercept)" "X" "C" "(Intercept)" ...
#>  @ mediator_predictors: chr [1:3] "(Intercept)" "X" "C"
#>  @ outcome_predictors : chr [1:4] "(Intercept)" "M" "X" "C"
#>  @ a_path             : num 0.367
#>  @ b_path             : num 0.343
#>  @ c_prime            : num 0.173
#>  @ vcov               : num [1:7, 1:7] 0.009243 -0.000939 0.001007 0 0 ...
#>  @ sigma_m            : num 0.951
#>  @ sigma_y            : num 1.05
#>  @ data               :'data.frame': 100 obs. of  3 variables:
#>  .. $ M: num  1.705 1.274 0.44 0.474 -0.635 ...
#>  .. $ X: num  -0.5605 -0.2302 1.5587 0.0705 0.1293 ...
#>  .. $ C: num  -0.71 0.257 -0.247 -0.348 -0.952 ...
#>  .. - attr(*, "terms")=Classes 'terms', 'formula'  language M ~ X + C
#>  ..  .. ..- attr(*, "variables")= language list(M, X, C)
#>  ..  .. ..- attr(*, "factors")= int [1:3, 1:2] 0 1 0 0 0 1
#>  ..  .. .. ..- attr(*, "dimnames")=List of 2
#>  ..  .. .. .. ..$ : chr [1:3] "M" "X" "C"
#>  ..  .. .. .. ..$ : chr [1:2] "X" "C"
#>  ..  .. ..- attr(*, "term.labels")= chr [1:2] "X" "C"
#>  ..  .. ..- attr(*, "order")= int [1:2] 1 1
#>  ..  .. ..- attr(*, "intercept")= int 1
#>  ..  .. ..- attr(*, "response")= int 1
#>  ..  .. ..- attr(*, ".Environment")=<environment: 0x55f96925a3e0> 
#>  ..  .. ..- attr(*, "predvars")= language list(M, X, C)
#>  ..  .. ..- attr(*, "dataClasses")= Named chr [1:3] "numeric" "numeric" "numeric"
#>  ..  .. .. ..- attr(*, "names")= chr [1:3] "M" "X" "C"
#>  @ n_obs              : int 100
#>  @ source_package     : chr "stats"
#>  @ converged          : logi TRUE
#>  @ treatment          : chr "X"
#>  @ mediator           : chr "M"
#>  @ outcome            : chr "Y"
# }
```
