# Integration with the mediation Package

The `probmed` package can extract mediation structures directly from
objects created by the popular `mediation` package. This allows you to
easily compute $P_{med}$ for models you have already fitted using
[`mediation::mediate()`](https://rdrr.io/pkg/mediation/man/mediate.html).

## Prerequisites

Ensure you have both `probmed` and `mediation` installed:

``` r
install.packages("mediation")
# install.packages("probmed") # Not on CRAN yet
```

Load the packages:

``` r
library(probmed)
library(mediation)
#> Loading required package: MASS
#> Loading required package: Matrix
#> Loading required package: mvtnorm
#> Loading required package: sandwich
#> mediation: Causal Mediation Analysis
#> Version: 4.5.1
```

## Basic Usage

The workflow is simple: 1. Fit your mediator and outcome models (using
`lm`, `glm`, etc.). 2. Run
[`mediation::mediate()`](https://rdrr.io/pkg/mediation/man/mediate.html).
3. Pass the resulting object to
[`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md).
4. Compute
[`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md).

### Example: Linear Mediation

``` r
# Simulate data
set.seed(123)
n <- 100
data <- data.frame(X = rnorm(n), C = rnorm(n))
data$M <- 0.5 * data$X + 0.3 * data$C + rnorm(n)
data$Y <- 0.4 * data$M + 0.2 * data$X + 0.1 * data$C + rnorm(n)

# 1. Fit models
model_m <- lm(M ~ X + C, data = data)
model_y <- lm(Y ~ M + X + C, data = data)

# 2. Run mediate
# We use boot=FALSE here for speed, but you can use bootstrap in mediate()
med_out <- mediate(model_m, model_y, treat = "X", mediator = "M", boot = FALSE)
summary(med_out)
#> 
#> Causal Mediation Analysis 
#> 
#> Quasi-Bayesian Confidence Intervals
#> 
#>                 Estimate 95% CI Lower 95% CI Upper p-value    
#> ACME            0.124547     0.034003     0.242980  <2e-16 ***
#> ADE             0.173342    -0.054517     0.398832   0.152    
#> Total Effect    0.297889     0.069188     0.528637   0.016 *  
#> Prop. Mediated  0.402200     0.100021     1.401127   0.016 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Sample Size Used: 100 
#> 
#> 
#> Simulations: 1000

# 3. Extract mediation structure
extract <- extract_mediation(med_out)
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
#>  ..  .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
#>  ..  .. ..- attr(*, "predvars")= language list(M, X, C)
#>  ..  .. ..- attr(*, "dataClasses")= Named chr [1:3] "numeric" "numeric" "numeric"
#>  ..  .. .. ..- attr(*, "names")= chr [1:3] "M" "X" "C"
#>  @ n_obs              : int 100
#>  @ source_package     : chr "stats"
#>  @ converged          : logi TRUE
#>  @ treatment          : chr "X"
#>  @ mediator           : chr "M"
#>  @ outcome            : chr "Y"

# 4. Compute P_med
pmed(extract, method = "parametric_bootstrap", n_boot = 200)
#> <probmed::PmedResult>
#>  @ estimate         : num 0.456
#>  @ ci_lower         : Named num 0.387
#>  .. - attr(*, "names")= chr "2.5%"
#>  @ ci_upper         : Named num 0.529
#>  .. - attr(*, "names")= chr "97.5%"
#>  @ ci_level         : num 0.95
#>  @ method           : chr "parametric_bootstrap"
#>  @ n_boot           : int 200
#>  @ boot_estimates   : num [1:200] 0.469 0.434 0.454 0.469 0.457 ...
#>  @ ie_estimate      : num 0.131
#>  @ ie_ci_lower      : Named num 0.0353
#>  .. - attr(*, "names")= chr "2.5%"
#>  @ ie_ci_upper      : Named num 0.267
#>  .. - attr(*, "names")= chr "97.5%"
#>  @ ie_boot_estimates: num [1:200] 0.073 0.08 0.14 0.198 0.267 ...
#>  @ x_ref            : num 0
#>  @ x_value          : num 1
#>  @ source_extract   : <probmed::MediationExtract>
#>  .. @ estimates          : Named num [1:7] 0.1351 0.3668 0.3238 -0.0193 0.3426 ...
#>  .. .. - attr(*, "names")= chr [1:7] "(Intercept)" "X" "C" "(Intercept)" ...
#>  .. @ mediator_predictors: chr [1:3] "(Intercept)" "X" "C"
#>  .. @ outcome_predictors : chr [1:4] "(Intercept)" "M" "X" "C"
#>  .. @ a_path             : num 0.367
#>  .. @ b_path             : num 0.343
#>  .. @ c_prime            : num 0.173
#>  .. @ vcov               : num [1:7, 1:7] 0.009243 -0.000939 0.001007 0 0 ...
#>  .. @ sigma_m            : num 0.951
#>  .. @ sigma_y            : num 1.05
#>  .. @ data               :'data.frame':  100 obs. of  3 variables:
#>  .. .. $ M: num  1.705 1.274 0.44 0.474 -0.635 ...
#>  .. .. $ X: num  -0.5605 -0.2302 1.5587 0.0705 0.1293 ...
#>  .. .. $ C: num  -0.71 0.257 -0.247 -0.348 -0.952 ...
#>  .. .. - attr(*, "terms")=Classes 'terms', 'formula'  language M ~ X + C
#>  .. ..  .. ..- attr(*, "variables")= language list(M, X, C)
#>  .. ..  .. ..- attr(*, "factors")= int [1:3, 1:2] 0 1 0 0 0 1
#>  .. ..  .. .. ..- attr(*, "dimnames")=List of 2
#>  .. ..  .. .. .. ..$ : chr [1:3] "M" "X" "C"
#>  .. ..  .. .. .. ..$ : chr [1:2] "X" "C"
#>  .. ..  .. ..- attr(*, "term.labels")= chr [1:2] "X" "C"
#>  .. ..  .. ..- attr(*, "order")= int [1:2] 1 1
#>  .. ..  .. ..- attr(*, "intercept")= int 1
#>  .. ..  .. ..- attr(*, "response")= int 1
#>  .. ..  .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
#>  .. ..  .. ..- attr(*, "predvars")= language list(M, X, C)
#>  .. ..  .. ..- attr(*, "dataClasses")= Named chr [1:3] "numeric" "numeric" "numeric"
#>  .. ..  .. .. ..- attr(*, "names")= chr [1:3] "M" "X" "C"
#>  .. @ n_obs              : int 100
#>  .. @ source_package     : chr "stats"
#>  .. @ converged          : logi TRUE
#>  .. @ treatment          : chr "X"
#>  .. @ mediator           : chr "M"
#>  .. @ outcome            : chr "Y"
#>  @ converged        : logi TRUE
#>  @ call             : NULL
```

## Why use probmed with mediation?

While the `mediation` package provides excellent tools for estimating
natural direct and indirect effects (ADE and ACME), `probmed` adds the
**Probabilistic Index of Mediation ($P_{med}$)**.

$P_{med}$ is a scale-free effect size that represents the probability
that the outcome for a treated individual is higher than for a control
individual, specifically due to the indirect path. This complements the
raw effect estimates provided by `mediation`.

## Supported Models

Since `probmed` extracts the underlying `lm` or `glm` objects from the
`mediate` object, it supports any model type that `probmed` supports
(currently Linear Models and Generalized Linear Models).

### Example: Binary Outcome (GLM)

``` r
# Simulate binary outcome
data$Y_prob <- plogis(0.4 * data$M + 0.2 * data$X)
data$Y_bin <- rbinom(n, 1, data$Y_prob)

# Fit models
model_m <- lm(M ~ X + C, data = data)
model_y <- glm(Y_bin ~ M + X + C, data = data, family = binomial)

# Run mediate
med_out_bin <- mediate(model_m, model_y, treat = "X", mediator = "M", boot = FALSE)

# Extract and compute P_med
extract_bin <- extract_mediation(med_out_bin)
pmed(extract_bin, method = "parametric_bootstrap", n_boot = 200)
#> <probmed::PmedResult>
#>  @ estimate         : num 0.065
#>  @ ci_lower         : Named num 0
#>  .. - attr(*, "names")= chr "2.5%"
#>  @ ci_upper         : Named num 1
#>  .. - attr(*, "names")= chr "97.5%"
#>  @ ci_level         : num 0.95
#>  @ method           : chr "parametric_bootstrap"
#>  @ n_boot           : int 200
#>  @ boot_estimates   : num [1:200] 0 0 0 0 0 0 0 0 0 0 ...
#>  @ ie_estimate      : num 0.0782
#>  @ ie_ci_lower      : Named num -0.0643
#>  .. - attr(*, "names")= chr "2.5%"
#>  @ ie_ci_upper      : Named num 0.282
#>  .. - attr(*, "names")= chr "97.5%"
#>  @ ie_boot_estimates: num [1:200] 0.0659 0.0928 0.0873 0.0323 0.0554 ...
#>  @ x_ref            : num 0
#>  @ x_value          : num 1
#>  @ source_extract   : <probmed::MediationExtract>
#>  .. @ estimates          : Named num [1:7] 0.135 0.367 0.324 -0.11 0.235 ...
#>  .. .. - attr(*, "names")= chr [1:7] "(Intercept)" "X" "C" "(Intercept)" ...
#>  .. @ mediator_predictors: chr [1:3] "(Intercept)" "X" "C"
#>  .. @ outcome_predictors : chr [1:4] "(Intercept)" "M" "X" "C"
#>  .. @ a_path             : num 0.367
#>  .. @ b_path             : num 0.235
#>  .. @ c_prime            : num 0.369
#>  .. @ vcov               : num [1:7, 1:7] 0.009243 -0.000939 0.001007 0 0 ...
#>  .. @ sigma_m            : num 0.951
#>  .. @ sigma_y            : num NA
#>  .. @ data               :'data.frame':  100 obs. of  3 variables:
#>  .. .. $ M: num  1.705 1.274 0.44 0.474 -0.635 ...
#>  .. .. $ X: num  -0.5605 -0.2302 1.5587 0.0705 0.1293 ...
#>  .. .. $ C: num  -0.71 0.257 -0.247 -0.348 -0.952 ...
#>  .. .. - attr(*, "terms")=Classes 'terms', 'formula'  language M ~ X + C
#>  .. ..  .. ..- attr(*, "variables")= language list(M, X, C)
#>  .. ..  .. ..- attr(*, "factors")= int [1:3, 1:2] 0 1 0 0 0 1
#>  .. ..  .. .. ..- attr(*, "dimnames")=List of 2
#>  .. ..  .. .. .. ..$ : chr [1:3] "M" "X" "C"
#>  .. ..  .. .. .. ..$ : chr [1:2] "X" "C"
#>  .. ..  .. ..- attr(*, "term.labels")= chr [1:2] "X" "C"
#>  .. ..  .. ..- attr(*, "order")= int [1:2] 1 1
#>  .. ..  .. ..- attr(*, "intercept")= int 1
#>  .. ..  .. ..- attr(*, "response")= int 1
#>  .. ..  .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
#>  .. ..  .. ..- attr(*, "predvars")= language list(M, X, C)
#>  .. ..  .. ..- attr(*, "dataClasses")= Named chr [1:3] "numeric" "numeric" "numeric"
#>  .. ..  .. .. ..- attr(*, "names")= chr [1:3] "M" "X" "C"
#>  .. @ n_obs              : int 100
#>  .. @ source_package     : chr "stats"
#>  .. @ converged          : logi TRUE
#>  .. @ treatment          : chr "X"
#>  .. @ mediator           : chr "M"
#>  .. @ outcome            : chr "Y_bin"
#>  @ converged        : logi TRUE
#>  @ call             : NULL
```
