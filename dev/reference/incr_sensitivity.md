# M–Y Confounding Sensitivity for the Incremental Mediated Elasticity

Per-`delta` M–Y unmeasured-confounding sensitivity for an
[IncrPmedResult](https://data-wise.github.io/probmed/dev/reference/IncrPmedResult.md).
For each tilt factor `delta`, reports the additive numerator bias that
zeroes the incremental mediated share (`tipping`) and, optionally, the
bias that drives `P_med^delta` to a user `threshold`. Reuses
[`pmed_sensitivity()`](https://data-wise.github.io/probmed/dev/reference/pmed_sensitivity.md);
the additive-offset model is exact for the Paper 2 (mean-based) share.

## Usage

``` r
incr_sensitivity(object, threshold = 0, ...)
```

## Arguments

- object:

  An
  [IncrPmedResult](https://data-wise.github.io/probmed/dev/reference/IncrPmedResult.md)
  from
  [`incr_pmed()`](https://data-wise.github.io/probmed/dev/reference/incr_pmed.md).

- threshold:

  Numeric or `NA`: a `P_med` value whose tipping bias to report per
  `delta` (default `0` = "bias that makes the share vanish").

- ...:

  Unused.

## Value

A data frame: one row per `delta` with columns `delta`, `med`, `tot`,
`Pmed`, `tipping`, `tipping_threshold`.

## Examples

``` r
set.seed(1)
n <- 800; C <- rnorm(n)
A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
M <- 0.6 * A + 0.4 * C + rnorm(n)
Y <- 0.5 * A + 0.7 * M + 0.3 * C + rnorm(n)
fit <- incr_pmed(data.frame(A, M, Y, C), deltas = c(0.5, 1, 2))
incr_sensitivity(fit, threshold = 0)
#>   delta        med       tot      Pmed     tipping tipping_threshold
#> 1   0.5 0.16221922 0.3872558 0.4188943 -0.16221922       -0.16221922
#> 2   1.0 0.09485621 0.2289208 0.4143626 -0.09485621       -0.09485621
#> 3   2.0 0.04693414 0.1132813 0.4143150 -0.04693414       -0.04693414
```
