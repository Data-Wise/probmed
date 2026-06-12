# Maximized log-likelihood under the constraint `a*b` = ie0.

Pins b = ie0 / a; Vy is free (its MLE is SSE_Y / n). Optimizes over (a,
log Vm); all linear nuisance terms are profiled by OLS.

## Usage

``` r
.mbco_ll_constrained_ie(prep, ie0)
```

## Value

Maximized joint log-likelihood, or NA_real\_ if the fit fails.
