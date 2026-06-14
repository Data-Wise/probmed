# Parametric bootstrap for parallel joint P_med.

Draws the structural coefficient block `(a_1..a_k, b_1..b_k)` from
`N(estimates, vcov)` (residual SDs held fixed, as in the single-mediator
parametric bootstrap) and evaluates the closed-form kernel per draw.

## Usage

``` r
.pmed_parallel_parametric_boot(
  extract,
  x_ref,
  x_value,
  n_boot,
  ci_level,
  seed,
  ...
)
```
