# Paired contrast of the proportion mediated between two scales

Contrasts `g` at two aggregation scales estimated from the **same**
clusters by differencing their influence functions cluster-by-cluster.
Because the shared sampling noise cancels, the paired standard error is
materially smaller than the naive one that treats the two scales as
independent.

## Usage

``` r
rg_flow_contrast(fit1, fit2, ci_level = 0.95)
```

## Arguments

- fit1, fit2:

  [RgFlowResult](https://data-wise.github.io/probmed/dev/reference/RgFlowResult.md)
  objects estimated on the same clusters. The contrast is
  `fit2@g - fit1@g`. Clusters present in only one of the two are
  dropped.

- ci_level:

  Numeric: confidence level.

## Value

List with `delta`, `se_paired`, `se_naive`, `gain`
(`se_naive / se_paired`), `z`, `p_value`, `ci`, and `n_clusters`.

## Details

The gain is largest for two nearby scales, where the influence functions
are most similar cluster-by-cluster; the endpoint contrast (unit scale
versus fully aggregated) is the worst case for pairing.

## Examples

``` r
set.seed(1)
d <- data.frame(
  A = rnorm(400), M = rnorm(400), Y = rnorm(400),
  sch = rep(1:40, each = 10)
)
f1 <- g_hat_eif(rg_coarse_grain(d, 1, cluster = "sch")$Z,
                rg_coarse_grain(d, 1, cluster = "sch")$cluster, m = 1)
cg2 <- rg_coarse_grain(d, 2, cluster = "sch")
f2 <- g_hat_eif(cg2$Z, cg2$cluster, m = 2)
rg_flow_contrast(f1, f2)$se_paired
#> [1] 0.5916062
```
