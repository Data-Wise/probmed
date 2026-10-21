# Coarse-grain multilevel data to aggregation scale `m`

Applies the coarse-graining operator `T_m`: within each cluster, units
are randomly partitioned into blocks of size `m` and each block is
replaced by its mean. The cluster structure is preserved, so every scale
keeps the same clusters and two scales remain comparable
cluster-by-cluster.

## Usage

``` r
rg_coarse_grain(
  data,
  m,
  treatment = "A",
  mediator = "M",
  outcome = "Y",
  cluster = "cluster"
)
```

## Arguments

- data:

  Data frame with the treatment, mediator, outcome, and cluster columns.

- m:

  Positive integer block size. `m = 1` returns the unit scale unchanged.

- treatment, mediator, outcome:

  Character: column names of `A`, `M`, `Y`.

- cluster:

  Character: column name of the cluster identifier.

## Value

List with `Z` (numeric matrix of coarse-grained units, columns `A`, `M`,
`Y`) and `cluster` (their cluster labels).

## Details

Blocks are complete by construction: within a cluster of size `n_s`, the
final incomplete block is discarded, so the number of coarse-grained
units is `sum(floor(n_s / m))`, not `n / m`. Clusters smaller than `m`
drop out entirely at that scale. The composition property
`T_{m2}(T_{m1}(.)) = T_{m1 * m2}(.)` holds exactly on the units a
partition retains.

This function calls [`sample.int()`](https://rdrr.io/r/base/sample.html)
and is therefore random. It never sets a seed; seed the caller (or wrap
in
[`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html))
for reproducibility.

## Examples

``` r
set.seed(1)
d <- data.frame(
  A = rnorm(200), M = rnorm(200), Y = rnorm(200),
  sch = rep(1:20, each = 10)
)
cg <- rg_coarse_grain(d, m = 5, cluster = "sch")
nrow(cg$Z)  # 40 = 20 clusters * floor(10 / 5)
#> [1] 40
```
