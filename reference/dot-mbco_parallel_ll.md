# Joint Gaussian log-likelihood for the parallel free model at given params.

Each mediator residual `M[,j] - a_j X` is OLS-profiled on
`Dm_list[[j]]`; the outcome residual `Y - M b` is OLS-profiled on `Dy`.
Returns the summed maximized log-likelihood for the supplied
`(a, Vm, b, Vy)`.

## Usage

``` r
.mbco_parallel_ll(prep, a, Vm, b, Vy)
```
