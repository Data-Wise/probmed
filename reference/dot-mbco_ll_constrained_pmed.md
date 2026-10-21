# Maximized log-likelihood under the constraint P_med = pnorm(qstar).

Profiles out intercepts, the direct effect, and covariate coefficients
by OLS, with both compatible coefficient signs explored. Positive Vm and
Vy/b^2 coordinates enforce the original signed constraint; b is profiled
by its quadratic stationary equation. Two numerical coordinates remain.

## Usage

``` r
.mbco_ll_constrained_pmed(prep, qstar, a_sign)
```

## Value

Maximized joint log-likelihood, or NA_real\_ if the fit fails.
