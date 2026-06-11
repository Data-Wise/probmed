# Maximized log-likelihood under the constraint P_med = pnorm(qstar).

Profiles out intercepts, the direct effect, and covariate coefficients
by OLS, so only (a, log Vm, log\|b\|) are optimized numerically (b is
sign-pinned to the target branch; Vy is fixed by the P_med constraint).
Nelder-Mead is used because the Vy \<= 0 infeasible region is a hard
penalty wall that defeats gradient methods.

## Usage

``` r
.mbco_ll_constrained_pmed(prep, qstar, a_sign)
```

## Value

Maximized joint log-likelihood, or NA_real\_ if the fit fails.
