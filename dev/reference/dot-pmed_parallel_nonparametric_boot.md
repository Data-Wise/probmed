# Nonparametric bootstrap for parallel joint P_med.

Resamples rows, refits the k mediator models and the outcome model, and
evaluates the closed-form kernel per resample.

## Usage

``` r
.pmed_parallel_nonparametric_boot(
  extract,
  x_ref,
  x_value,
  n_boot,
  ci_level,
  seed,
  ...
)
```
