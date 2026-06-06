# Compute P_med: Probabilistic Effect Size for Mediation

Compute \\P\_{med}\\, a scale-free probabilistic effect size for
mediation analysis, along with the traditional Indirect Effect (\\a
\times b\\). Provides point estimates and bootstrap confidence
intervals.

\\P\_{med}\\ represents \\P(Y(x, M(x)) \> Y(x, M(x^\*))) + \tfrac12
P(Y(x, M(x)) = Y(x, M(x^\*)))\\: the probability that the outcome under
treatment with the mediator at its treated level exceeds the outcome
under treatment with the mediator at its control level (a one-half
correction is added for ties). Both potential outcomes hold treatment at
\\x\\, so the direct effect cancels and \\P\_{med}\\ reflects mediation
only; values \\\> 0.5\\ indicate positive mediation.

## Usage

``` r
pmed(object, ...)
```

## Arguments

- object:

  Either:

  - A `formula` for the outcome model (most common)

  - A
    [`medfit::MediationData`](https://data-wise.github.io/medfit/reference/MediationData.html)
    object from
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
#> 
#> P_med: Probability of Mediated Shift
#> ====================================
#> 
#> Estimate: 0.551 
#> 95% CI: [ 0.529 ,  0.578 ]
#> 
#> Indirect Effect (Product of Coefficients):
#> Estimate: 0.204 
#> 95% CI: [ 0.116 ,  0.317 ]
#> 
#> Inference: parametric_bootstrap 
#> Bootstrap samples: 500 
#> 
#> Treatment contrast: X = 1 vs. X* = 0 
#> 
#> Interpretation:
#>   P(Y_{X*, M_X} > Y_{X, M_X}) = 0.551 
#> 

if (FALSE) { # \dontrun{
# View bootstrap distribution
summary(result)
plot(result)
} # }
```
