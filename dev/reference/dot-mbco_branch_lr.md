# Branch likelihood-ratio statistics for H0: a = 0 and H0: b = 0

Nested-model LR statistics with ML variances and the linear nuisance
terms profiled, matching
[`.mbco_ll_null()`](https://data-wise.github.io/probmed/dev/reference/dot-mbco_ll_null.md).
When both are below the interval's chi-square(1) cutoff, neither path is
individually significant: the data sit near the double null a = b = 0,
where the profile LR at p\* = 0.5 is the minimum of these two statistics
and the chi-square(1) interval is conservative (coverage about 1 -
alpha^2, 0.9975 at the 95% level).

## Usage

``` r
.mbco_branch_lr(prep)
```

## Value

Named numeric vector `c(a0 = , b0 = )`.
