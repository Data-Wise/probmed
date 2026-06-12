# Run Nelder-Mead from several starts, returning the best maximized log-likelihood (`-min value`), or NA if every start fails / is penalized.

The first start (the analytic warm start) gets the full iteration
budget; the remaining (perturbed) starts are capped low – they exist
only to escape a nearby false basin, so a start that does not converge
quickly is abandoned cheaply rather than burning the full budget against
the penalty wall.

## Usage

``` r
.mbco_parallel_optim_best(nll, starts, maxit = 1500L, maxit_extra = 400L)
```
