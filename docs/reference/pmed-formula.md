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
#> <probmed::PmedResult>
#>  @ estimate         : num 0.391
#>  @ ci_lower         : num NA
#>  @ ci_upper         : num NA
#>  @ ci_level         : num NA
#>  @ method           : chr "plugin"
#>  @ n_boot           : int NA
#>  @ boot_estimates   : num(0) 
#>  @ ie_estimate      : num 0.148
#>  @ ie_ci_lower      : num NA
#>  @ ie_ci_upper      : num NA
#>  @ ie_boot_estimates: num(0) 
#>  @ x_ref            : num 0
#>  @ x_value          : num 1
#>  @ source_extract   : <probmed::MediationExtract>
#>  .. @ estimates          : Named num [1:7] 0.06753 0.43341 0.31191 -0.00966 0.20097 ...
#>  .. .. - attr(*, "names")= chr [1:7] "(Intercept)" "X" "C" "(Intercept)" ...
#>  .. @ mediator_predictors: chr [1:3] "(Intercept)" "X" "C"
#>  .. @ outcome_predictors : chr [1:4] "(Intercept)" "X" "M" "C"
#>  .. @ a_path             : num 0.433
#>  .. @ b_path             : num 0.343
#>  .. @ c_prime            : num 0.201
#>  .. @ vcov               : num [1:7, 1:7] 0.002311 -0.000235 0.000252 0 0 ...
#>  .. @ sigma_m            : num 0.476
#>  .. @ sigma_y            : num 0.526
#>  .. @ data               :'data.frame':  100 obs. of  4 variables:
#>  .. .. $ X: num  -0.5605 -0.2302 1.5587 0.0705 0.1293 ...
#>  .. .. $ C: num  -0.71 0.257 -0.247 -0.348 -0.952 ...
#>  .. .. $ M: num  0.606 0.618 0.573 0.203 -0.428 ...
#>  .. .. $ Y: num  -0.3694 -0.1237 0.0222 -0.5006 -0.5543 ...
#>  .. @ n_obs              : int 100
#>  .. @ source_package     : chr "stats"
#>  .. @ converged          : logi TRUE
#>  .. @ treatment          : chr "X"
#>  .. @ mediator           : chr "M"
#>  .. @ outcome            : chr "Y"
#>  @ converged        : logi TRUE
#>  @ call             : language `method(pmed, new_S3_class("formula"))`(object = Y ~ X + M + C, formula_m = ..1,      data = ..2, treatment = "X"| __truncated__

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
#> <probmed::PmedResult>
#>  @ estimate         : num 0.397
#>  @ ci_lower         : Named num 0.322
#>  .. - attr(*, "names")= chr "2.5%"
#>  @ ci_upper         : Named num 0.462
#>  .. - attr(*, "names")= chr "97.5%"
#>  @ ci_level         : num 0.95
#>  @ method           : chr "parametric_bootstrap"
#>  @ n_boot           : int 200
#>  @ boot_estimates   : num [1:200] 0.362 0.418 0.424 0.321 0.416 ...
#>  @ ie_estimate      : num 0.146
#>  @ ie_ci_lower      : Named num 0.053
#>  .. - attr(*, "names")= chr "2.5%"
#>  @ ie_ci_upper      : Named num 0.242
#>  .. - attr(*, "names")= chr "97.5%"
#>  @ ie_boot_estimates: num [1:200] 0.0806 0.2042 0.1909 0.112 0.0958 ...
#>  @ x_ref            : num 0
#>  @ x_value          : num 1
#>  @ source_extract   : <probmed::MediationExtract>
#>  .. @ estimates          : Named num [1:7] 0.06753 0.43341 0.31191 -0.00966 0.20097 ...
#>  .. .. - attr(*, "names")= chr [1:7] "(Intercept)" "X" "C" "(Intercept)" ...
#>  .. @ mediator_predictors: chr [1:3] "(Intercept)" "X" "C"
#>  .. @ outcome_predictors : chr [1:4] "(Intercept)" "X" "M" "C"
#>  .. @ a_path             : num 0.433
#>  .. @ b_path             : num 0.343
#>  .. @ c_prime            : num 0.201
#>  .. @ vcov               : num [1:7, 1:7] 0.002311 -0.000235 0.000252 0 0 ...
#>  .. @ sigma_m            : num 0.476
#>  .. @ sigma_y            : num 0.526
#>  .. @ data               :'data.frame':  100 obs. of  4 variables:
#>  .. .. $ X: num  -0.5605 -0.2302 1.5587 0.0705 0.1293 ...
#>  .. .. $ C: num  -0.71 0.257 -0.247 -0.348 -0.952 ...
#>  .. .. $ M: num  0.606 0.618 0.573 0.203 -0.428 ...
#>  .. .. $ Y: num  -0.3694 -0.1237 0.0222 -0.5006 -0.5543 ...
#>  .. @ n_obs              : int 100
#>  .. @ source_package     : chr "stats"
#>  .. @ converged          : logi TRUE
#>  .. @ treatment          : chr "X"
#>  .. @ mediator           : chr "M"
#>  .. @ outcome            : chr "Y"
#>  @ converged        : logi TRUE
#>  @ call             : language `method(pmed, new_S3_class("formula"))`(object = Y ~ X + M + C, formula_m = ..1,      data = ..2, treatment = "X"| __truncated__ ...
# }
```
