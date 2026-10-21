# Doubly-Robust Wasserstein Proportion Mediated

Standalone cross-fit DR estimator of `P_med^W`. Consistent when EITHER
the outcome model OR the propensity model is misspecified, provided the
mediator model is correct (multiply robust; verified in simulation).

## Usage

``` r
pmedW_dr(d, covars = "C", K = 5L, misspec = "none", G = 300L, seed = 1L)
```

## Arguments

- d:

  A `data.frame` with columns `A`, `M`, `Y`, and the covariates in
  `covars`.

- covars:

  Character vector of covariate column names. Default `"C"`.

- K:

  Integer: cross-fitting folds (default 5).

- misspec:

  Character: `"none"` (both models correct), `"outcome"`, or
  `"propensity"` — for robustness testing only.

- G:

  Integer: number of `y`-grid points for the CDF estimator (default
  300).

- seed:

  Integer: RNG seed for fold assignment (default 1).

## Value

A named list with `NIE_W`, `NDE_W`, and `pmedW`.
