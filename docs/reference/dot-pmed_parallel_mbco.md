# MBCO confidence interval for parallel joint P_med (Gaussian).

Deterministic likelihood-ratio inversion (Tofighi & Kelley, 2020) for
the joint parallel estimand and for the total indirect effect
`sum(a*b)`, reusing the generic `.mbco_invert` engine. Gaussian outcome
and mediators only.

## Usage

``` r
.pmed_parallel_mbco(extract, x_ref, x_value, ci_level = 0.95, ...)
```
