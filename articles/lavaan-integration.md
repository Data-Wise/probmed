# Working with lavaan Models

The `probmed` package integrates seamlessly with `lavaan`, the popular R
package for structural equation modeling (SEM). This vignette
demonstrates how to extract mediation structures from `lavaan` models
and compute the probabilistic index of mediation ($P_{med}$).

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
#> This is lavaan 0.6-20
#> lavaan is FREE software! Please report any bugs.
```

## Simple Mediation Model

Let’s start with a simple mediation model where $X$ affects $M$, and
both $X$ and $M$ affect $Y$.

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
#> lavaan 0.6-20 ended normally after 1 iteration
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                         5
#> 
#>   Number of observations                           200
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 0
#> 
#> Parameter Estimates:
#> 
#>   Standard errors                             Standard
#>   Information                                 Expected
#>   Information saturated (h1) model          Structured
#> 
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   M ~                                                 
#>     X          (a)    0.471    0.075    6.307    0.000
#>   Y ~                                                 
#>     M          (b)    0.349    0.068    5.107    0.000
#>     X         (cp)    0.193    0.079    2.440    0.015
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .M                 0.986    0.099   10.000    0.000
#>    .Y                 0.923    0.092   10.000    0.000
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
#> <probmed::MediationExtract>
#>  @ estimates          : Named num [1:3] 0.471 0.349 0.193
#>  .. - attr(*, "names")= chr [1:3] "a" "b" "cp"
#>  @ mediator_predictors: chr "X"
#>  @ outcome_predictors : chr [1:2] "M" "X"
#>  @ a_path             : num 0.471
#>  @ b_path             : num 0.349
#>  @ c_prime            : num 0.193
#>  @ vcov               : num [1:3, 1:3] 5.57e-03 -5.07e-19 -2.24e-19 -6.51e-19 4.68e-03 ...
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:3] "a" "b" "cp"
#>  ..  ..$ : chr [1:3] "a" "b" "cp"
#>  @ sigma_m            : num 0.993
#>  @ sigma_y            : num 0.961
#>  @ data               :'data.frame': 200 obs. of  3 variables:
#>  .. $ M: num  1.919 1.197 0.514 0.578 -0.35 ...
#>  .. $ Y: num  0.582 -0.736 -0.117 0.217 0.557 ...
#>  .. $ X: num  -0.5605 -0.2302 1.5587 0.0705 0.1293 ...
#>  @ n_obs              : int 200
#>  @ source_package     : chr "lavaan"
#>  @ converged          : logi TRUE
#>  @ treatment          : chr "X"
#>  @ mediator           : chr "M"
#>  @ outcome            : chr "Y"
```

### 4. Compute $P_{med}$

Now you can compute $P_{med}$ using any of the supported methods
(plugin, parametric bootstrap, or non-parametric bootstrap).

``` r
# Plugin estimator (fast, point estimate only)
pmed(extract, method = "plugin")
#> <probmed::PmedResult>
#>  @ estimate         : num 0.441
#>  @ ci_lower         : num NA
#>  @ ci_upper         : num NA
#>  @ ci_level         : num NA
#>  @ method           : chr "plugin"
#>  @ n_boot           : int NA
#>  @ boot_estimates   : num(0) 
#>  @ ie_estimate      : num 0.164
#>  @ ie_ci_lower      : num NA
#>  @ ie_ci_upper      : num NA
#>  @ ie_boot_estimates: num(0) 
#>  @ x_ref            : num 0
#>  @ x_value          : num 1
#>  @ source_extract   : <probmed::MediationExtract>
#>  .. @ estimates          : Named num [1:3] 0.471 0.349 0.193
#>  .. .. - attr(*, "names")= chr [1:3] "a" "b" "cp"
#>  .. @ mediator_predictors: chr "X"
#>  .. @ outcome_predictors : chr [1:2] "M" "X"
#>  .. @ a_path             : num 0.471
#>  .. @ b_path             : num 0.349
#>  .. @ c_prime            : num 0.193
#>  .. @ vcov               : num [1:3, 1:3] 5.57e-03 -5.07e-19 -2.24e-19 -6.51e-19 4.68e-03 ...
#>  .. .. - attr(*, "dimnames")=List of 2
#>  .. ..  ..$ : chr [1:3] "a" "b" "cp"
#>  .. ..  ..$ : chr [1:3] "a" "b" "cp"
#>  .. @ sigma_m            : num 0.993
#>  .. @ sigma_y            : num 0.961
#>  .. @ data               :'data.frame':  200 obs. of  3 variables:
#>  .. .. $ M: num  1.919 1.197 0.514 0.578 -0.35 ...
#>  .. .. $ Y: num  0.582 -0.736 -0.117 0.217 0.557 ...
#>  .. .. $ X: num  -0.5605 -0.2302 1.5587 0.0705 0.1293 ...
#>  .. @ n_obs              : int 200
#>  .. @ source_package     : chr "lavaan"
#>  .. @ converged          : logi TRUE
#>  .. @ treatment          : chr "X"
#>  .. @ mediator           : chr "M"
#>  .. @ outcome            : chr "Y"
#>  @ converged        : logi TRUE
#>  @ call             : NULL

# Parametric bootstrap (with confidence intervals)
# Note: Increase n_boot for production use
pmed(extract, method = "parametric_bootstrap", n_boot = 200, seed = 456)
#> <probmed::PmedResult>
#>  @ estimate         : num 0.445
#>  @ ci_lower         : Named num 0.405
#>  .. - attr(*, "names")= chr "2.5%"
#>  @ ci_upper         : Named num 0.494
#>  .. - attr(*, "names")= chr "97.5%"
#>  @ ci_level         : num 0.95
#>  @ method           : chr "parametric_bootstrap"
#>  @ n_boot           : int 200
#>  @ boot_estimates   : num [1:200] 0.478 0.415 0.429 0.465 0.454 ...
#>  @ ie_estimate      : num 0.164
#>  @ ie_ci_lower      : Named num 0.102
#>  .. - attr(*, "names")= chr "2.5%"
#>  @ ie_ci_upper      : Named num 0.251
#>  .. - attr(*, "names")= chr "97.5%"
#>  @ ie_boot_estimates: num [1:200] 0.19 0.185 0.146 0.203 0.179 ...
#>  @ x_ref            : num 0
#>  @ x_value          : num 1
#>  @ source_extract   : <probmed::MediationExtract>
#>  .. @ estimates          : Named num [1:3] 0.471 0.349 0.193
#>  .. .. - attr(*, "names")= chr [1:3] "a" "b" "cp"
#>  .. @ mediator_predictors: chr "X"
#>  .. @ outcome_predictors : chr [1:2] "M" "X"
#>  .. @ a_path             : num 0.471
#>  .. @ b_path             : num 0.349
#>  .. @ c_prime            : num 0.193
#>  .. @ vcov               : num [1:3, 1:3] 5.57e-03 -5.07e-19 -2.24e-19 -6.51e-19 4.68e-03 ...
#>  .. .. - attr(*, "dimnames")=List of 2
#>  .. ..  ..$ : chr [1:3] "a" "b" "cp"
#>  .. ..  ..$ : chr [1:3] "a" "b" "cp"
#>  .. @ sigma_m            : num 0.993
#>  .. @ sigma_y            : num 0.961
#>  .. @ data               :'data.frame':  200 obs. of  3 variables:
#>  .. .. $ M: num  1.919 1.197 0.514 0.578 -0.35 ...
#>  .. .. $ Y: num  0.582 -0.736 -0.117 0.217 0.557 ...
#>  .. .. $ X: num  -0.5605 -0.2302 1.5587 0.0705 0.1293 ...
#>  .. @ n_obs              : int 200
#>  .. @ source_package     : chr "lavaan"
#>  .. @ converged          : logi TRUE
#>  .. @ treatment          : chr "X"
#>  .. @ mediator           : chr "M"
#>  .. @ outcome            : chr "Y"
#>  @ converged        : logi TRUE
#>  @ call             : NULL
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
#> <probmed::PmedResult>
#>  @ estimate         : num 0.443
#>  @ ci_lower         : num NA
#>  @ ci_upper         : num NA
#>  @ ci_level         : num NA
#>  @ method           : chr "plugin"
#>  @ n_boot           : int NA
#>  @ boot_estimates   : num(0) 
#>  @ ie_estimate      : num 0.162
#>  @ ie_ci_lower      : num NA
#>  @ ie_ci_upper      : num NA
#>  @ ie_boot_estimates: num(0) 
#>  @ x_ref            : num 0
#>  @ x_value          : num 1
#>  @ source_extract   : <probmed::MediationExtract>
#>  .. @ estimates          : Named num [1:3] 0.461 0.352 0.191
#>  .. .. - attr(*, "names")= chr [1:3] "a" "b" "cp"
#>  .. @ mediator_predictors: chr [1:3] "X" "C1" "C2"
#>  .. @ outcome_predictors : chr [1:4] "M" "X" "C1" "C2"
#>  .. @ a_path             : num 0.461
#>  .. @ b_path             : num 0.352
#>  .. @ c_prime            : num 0.191
#>  .. @ vcov               : num [1:3, 1:3] 5.60e-03 -3.71e-18 2.75e-18 -3.45e-18 4.67e-03 ...
#>  .. .. - attr(*, "dimnames")=List of 2
#>  .. ..  ..$ : chr [1:3] "a" "b" "cp"
#>  .. ..  ..$ : chr [1:3] "a" "b" "cp"
#>  .. @ sigma_m            : num 0.989
#>  .. @ sigma_y            : num 0.957
#>  .. @ data               :'data.frame':  200 obs. of  5 variables:
#>  .. .. $ M : num  1.919 1.197 0.514 0.578 -0.35 ...
#>  .. .. $ Y : num  0.582 -0.736 -0.117 0.217 0.557 ...
#>  .. .. $ X : num  -0.5605 -0.2302 1.5587 0.0705 0.1293 ...
#>  .. .. $ C1: num  -0.832 1.716 0.916 -0.46 -1.405 ...
#>  .. .. $ C2: num  -1.147 -0.441 -1.558 0.544 -0.233 ...
#>  .. @ n_obs              : int 200
#>  .. @ source_package     : chr "lavaan"
#>  .. @ converged          : logi TRUE
#>  .. @ treatment          : chr "X"
#>  .. @ mediator           : chr "M"
#>  .. @ outcome            : chr "Y"
#>  @ converged        : logi TRUE
#>  @ call             : NULL
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
#> Warning in `method(extract_mediation, lavaan)`(object = new("lavaan", version =
#> "0.6-20", : Variance-covariance matrix corresponds to unstandardized estimates.
#> Standard errors for standardized estimates may be incorrect if used directly.

print(extract_std)
#> <probmed::MediationExtract>
#>  @ estimates          : Named num [1:3] 0.407 0.353 0.169
#>  .. - attr(*, "names")= chr [1:3] "a" "b" "cp"
#>  @ mediator_predictors: chr "X"
#>  @ outcome_predictors : chr [1:2] "M" "X"
#>  @ a_path             : num 0.407
#>  @ b_path             : num 0.353
#>  @ c_prime            : num 0.169
#>  @ vcov               : num [1:3, 1:3] 5.57e-03 -5.07e-19 -2.24e-19 -6.51e-19 4.68e-03 ...
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:3] "a" "b" "cp"
#>  ..  ..$ : chr [1:3] "a" "b" "cp"
#>  @ sigma_m            : num 0.913
#>  @ sigma_y            : num 0.893
#>  @ data               :'data.frame': 200 obs. of  3 variables:
#>  .. $ M: num  1.919 1.197 0.514 0.578 -0.35 ...
#>  .. $ Y: num  0.582 -0.736 -0.117 0.217 0.557 ...
#>  .. $ X: num  -0.5605 -0.2302 1.5587 0.0705 0.1293 ...
#>  @ n_obs              : int 200
#>  @ source_package     : chr "lavaan"
#>  @ converged          : logi TRUE
#>  @ treatment          : chr "X"
#>  @ mediator           : chr "M"
#>  @ outcome            : chr "Y"
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

Ensure that the model actually contains the paths $M \sim X$,
$Y \sim M$, and $Y \sim X$. If any of these are missing (e.g., if you
fixed $c\prime$ to 0),
[`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md)
will fail because $P_{med}$ requires all three paths (though $c\prime$
can be 0, it must be in the model structure or handled as 0).

Currently, `probmed` expects the direct effect $Y \sim X$ to be present
in the model. If you omitted it (assuming full mediation), you can add
it and fix it to 0 in `lavaan`: `Y ~ 0*X`.
