# Multivariate Wasserstein Proportion Mediated

Computes `P_med^W` for multivariate outcomes (`d >= 1`) from corner
pseudo-samples via debiased Sinkhorn divergence. Intended for use with
pre-computed corner-law draws (e.g. from g-computation or DR). SE not
yet implemented.

## Usage

``` r
pmedW_md(NU11, NU10, NU00, eps = 0.1)
```

## Arguments

- NU11, NU10, NU00:

  Numeric matrices (rows = draws, columns = outcome dimensions).

- eps:

  Numeric: Sinkhorn regularization parameter (default 0.1).

## Value

A named list with `NIE_W`, `NDE_W`, and `pmedW`.
