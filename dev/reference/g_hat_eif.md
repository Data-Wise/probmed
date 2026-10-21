# Efficient-influence-function estimator of the scale-indexed proportion mediated

Estimates `g_l = NIE_l / TE_l` on coarse-grained data by plug-in, with a
**cluster-robust** standard error from the efficient influence function.
The plug-in is already efficient in the linear-Gaussian case, so no
cross-fitting is performed.

## Usage

``` r
g_hat_eif(Z, cluster, m = 1, ci_level = 0.95)
```

## Arguments

- Z:

  Numeric matrix with columns `A`, `M`, `Y` (coarse-grained units, as
  returned by
  [`rg_coarse_grain()`](https://data-wise.github.io/probmed/dev/reference/rg_coarse_grain.md)).

- cluster:

  Cluster labels for the rows of `Z`. The i.i.d. sampling unit is the
  cluster, not the row; passing unit-level labels would understate every
  standard error.

- m:

  Numeric: the aggregation scale `Z` was produced at (recorded, not
  recomputed).

- ci_level:

  Numeric: confidence level for the Wald interval.

## Value

An
[RgFlowResult](https://data-wise.github.io/probmed/dev/reference/RgFlowResult.md).

## Examples

``` r
set.seed(1)
d <- data.frame(
  A = rnorm(400), M = rnorm(400), Y = rnorm(400),
  sch = rep(1:40, each = 10)
)
cg <- rg_coarse_grain(d, m = 2, cluster = "sch")
fit <- g_hat_eif(cg$Z, cg$cluster, m = 2)
fit@g
#> [1] 0.03047159
```
