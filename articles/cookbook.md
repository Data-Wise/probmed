# Cookbook: common tasks

Short recipes for common tasks. Every recipe uses the simulated data
below (continuous treatment `X`, mediator `M`, outcome `Y`, covariate
`C`). Chunks are not evaluated when the site is built; the code was run
against version 0.4.0.

``` r

library(probmed)
set.seed(123)
n <- 300
d <- data.frame(X = rnorm(n), C = rnorm(n))
d$M <- 0.5 * d$X + 0.3 * d$C + rnorm(n)
d$Y <- 0.4 * d$M + 0.2 * d$X + 0.2 * d$C + rnorm(n)
```

## Point estimate only

``` r

fit <- pmed(Y ~ X + M + C, formula_m = M ~ X + C, data = d,
            treatment = "X", mediator = "M", method = "plugin", seed = 1)
fit@estimate
```

`"plugin"` returns the point estimate with no interval. Use it for quick
checks.

## Interval by method

All methods share one call; only `method` changes.

``` r

args <- list(Y ~ X + M + C, formula_m = M ~ X + C, data = d,
             treatment = "X", mediator = "M", seed = 1)
pb   <- do.call(pmed, c(args, method = "parametric_bootstrap", n_boot = 1000))
npb  <- do.call(pmed, c(args, method = "nonparametric_bootstrap", n_boot = 1000))
mbco <- do.call(pmed, c(args, method = "mbco"))
```

See [Choosing a confidence
interval](https://data-wise.github.io/probmed/articles/confidence-intervals.md)
for when to use which.

## Read the result

``` r

mbco@estimate; mbco@ci_lower; mbco@ci_upper   # P_med and its interval
mbco@ie_estimate                               # indirect effect a*b
mbco@diagnostics$pmed                          # interval status, near_double_null flag
```

For `method = "mbco"`, check that the printed interval status reads
`ok`. A missing bound means a verified likelihood-ratio endpoint was not
found; the method reports that instead of a guess.

## Binary mediator or outcome with covariates

Use the model-based bootstrap, which refits both equations on simulated
data.

``` r

pmed(Y_bin ~ X + M + C, formula_m = M ~ X + C, data = d_bin,
     treatment = "X", mediator = "M",
     method = "model_bootstrap", n_boot = 500, seed = 1)
```

A rare outcome (about 5%) at small N often returns
`bootstrap_incomplete` instead of an interval. That is the method
declining to report, not an error. Poisson, ordinal and other families,
interactions and transformed treatment or mediator terms are rejected.

## From a fitted lavaan or mediation object

[`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md)
returns a `MediationData` object, and
[`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md) has a
method for it. See [lavaan
integration](https://data-wise.github.io/probmed/articles/lavaan-integration.md)
and [mediation
integration](https://data-wise.github.io/probmed/articles/mediation-integration.md).

## Two mediators

For parallel mediators use the joint P_med; see [Parallel
mediation](https://data-wise.github.io/probmed/articles/parallel-mediation.md).

## Tipping point for unmeasured M-Y confounding (mean-based shares)

``` r

pmed_sensitivity(indirect = 0.30, total = 0.50, threshold = 0, evalue = TRUE)
```

More in [Sensitivity
analysis](https://data-wise.github.io/probmed/articles/sensitivity-analysis.md).

## Weak identification check

``` r

ward_residual(data.frame(A = as.integer(d$X > 0), M = d$M, Y = d$Y, C = d$C))
```

[`ward_residual()`](https://data-wise.github.io/probmed/reference/ward_residual.md)
expects a binary 0/1 treatment `A` and covariates named in `covars`
(default `"C"`). It reports gates for weak identification; see [Gauge
residual](https://data-wise.github.io/probmed/articles/gauge-residual.md).

## Reproducibility

Pass `seed` to every call. Plug-in calls honor it, MBCO is
deterministic, and bootstrap methods draw from `seed`.
