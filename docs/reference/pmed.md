# Compute P_med: Probabilistic Effect Size for Mediation

Compute \\P\_{med}\\, a scale-free probabilistic effect size for
mediation analysis, along with the traditional Indirect Effect (\\a
\times b\\). Provides point estimates and bootstrap confidence
intervals.

\\P\_{med}\\ represents \\P(Y\_{x^\*, M_x} \> Y\_{x, M_x})\\, the
probability that the counterfactual outcome under control with mediator
at treated level exceeds the outcome under treatment.

## Usage

``` r
pmed(object, ...)
```

## Arguments

- object:

  Either:

  - A `formula` for the outcome model (most common)

  - A `MediationExtract` object from
    [`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md)

- ...:

  Additional arguments passed to methods (see method documentation)

## Value

`PmedResult` object containing:

- `estimate`: P_med point estimate

- `ci_lower`, `ci_upper`: Confidence interval bounds

- `ie_estimate`: Indirect Effect point estimate

- `ie_ci_lower`, `ie_ci_upper`: IE confidence interval

- `boot_estimates`: Bootstrap distribution (if applicable)

- `method`: Inference method used

## Methods

Available inference methods via `method` argument:

- `"parametric_bootstrap"` (default): Fast, assumes normality

- `"nonparametric_bootstrap"`: Robust to assumptions, slower

- `"plugin"`: Point estimate only, no CI

## See also

[`extract_mediation`](https://data-wise.github.io/probmed/reference/extract_mediation.md)
for extracting from fitted models

## Examples

``` r
# Basic example with formula interface
set.seed(123)
n <- 200
data <- data.frame(X = rnorm(n), C = rnorm(n))
data$M <- 0.5 * data$X + 0.3 * data$C + rnorm(n)
data$Y <- 0.4 * data$M + 0.2 * data$X + 0.2 * data$C + rnorm(n)

# Compute with parametric bootstrap
result <- pmed(
  Y ~ X + M + C,
  formula_m = M ~ X + C,
  data = data,
  treatment = "X",
  mediator = "M",
  method = "parametric_bootstrap",
  n_boot = 500
)
print(result)
#> <probmed::PmedResult>
#>  @ estimate         : num 0.463
#>  @ ci_lower         : Named num 0.42
#>  .. - attr(*, "names")= chr "2.5%"
#>  @ ci_upper         : Named num 0.511
#>  .. - attr(*, "names")= chr "97.5%"
#>  @ ci_level         : num 0.95
#>  @ method           : chr "parametric_bootstrap"
#>  @ n_boot           : int 500
#>  @ boot_estimates   : num [1:500] 0.479 0.452 0.5 0.485 0.457 ...
#>  @ ie_estimate      : num 0.204
#>  @ ie_ci_lower      : Named num 0.122
#>  .. - attr(*, "names")= chr "2.5%"
#>  @ ie_ci_upper      : Named num 0.315
#>  .. - attr(*, "names")= chr "97.5%"
#>  @ ie_boot_estimates: num [1:500] 0.159 0.158 0.213 0.289 0.276 ...
#>  @ x_ref            : num 0
#>  @ x_value          : num 1
#>  @ source_extract   : <probmed::MediationExtract>
#>  .. @ estimates          : Named num [1:7] 0.0336 0.4676 0.2493 -0.027 0.1347 ...
#>  .. .. - attr(*, "names")= chr [1:7] "(Intercept)" "X" "C" "(Intercept)" ...
#>  .. @ mediator_predictors: chr [1:3] "(Intercept)" "X" "C"
#>  .. @ outcome_predictors : chr [1:4] "(Intercept)" "X" "M" "C"
#>  .. @ a_path             : num 0.468
#>  .. @ b_path             : num 0.437
#>  .. @ c_prime            : num 0.135
#>  .. @ vcov               : num [1:7, 1:7] 4.69e-03 3.95e-05 -1.99e-04 0.00 0.00 ...
#>  .. @ sigma_m            : num 0.968
#>  .. @ sigma_y            : num 1.03
#>  .. @ data               :'data.frame':  200 obs. of  4 variables:
#>  .. .. $ X: num  -0.5605 -0.2302 1.5587 0.0705 0.1293 ...
#>  .. .. $ C: num  2.199 1.312 -0.265 0.543 -0.414 ...
#>  .. .. $ M: num  0.3058 -0.89 0.0651 0.1694 0.611 ...
#>  .. .. $ Y: num  1.524 -0.167 0.251 -1.326 0.978 ...
#>  .. @ n_obs              : int 200
#>  .. @ source_package     : chr "stats"
#>  .. @ converged          : logi TRUE
#>  .. @ treatment          : chr "X"
#>  .. @ mediator           : chr "M"
#>  .. @ outcome            : chr "Y"
#>  @ converged        : logi TRUE
#>  @ call             : language `method(pmed, new_S3_class("formula"))`(object = Y ~ X + M + C, formula_m = ..1,      data = ..2, treatment = "X"| __truncated__ ...

if (FALSE) { # \dontrun{
# View bootstrap distribution
summary(result)
plot(result)
} # }
```
