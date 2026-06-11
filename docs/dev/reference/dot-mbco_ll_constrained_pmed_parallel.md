# Maximized log-likelihood under the constraint joint P_med = pnorm(qstar).

The constraint fixes `Vy = (delta*S)^2 / (2 qstar^2) - sum(b^2 Vm)` with
`S = sum(a*b)`, and requires `sign(delta*S) = sign(qstar)`. Optimizes
over `(a_1..a_k, log Vm_1..k, b_1..k)` (3k coords) with intercepts /
direct effect / covariates OLS-profiled. Nelder-Mead (infeasible Vy is a
penalty wall).

## Usage

``` r
.mbco_ll_constrained_pmed_parallel(prep, qstar)
```
