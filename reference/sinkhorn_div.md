# Debiased Sinkhorn Divergence

Debiased Sinkhorn divergence `S_eps(X, Y)` (Genevay et al.), which
approximates `W_2^2(X, Y)` as `eps -> 0`. Use for multivariate (`d > 1`)
Wasserstein P_med^W via
[`pmedW_md()`](https://data-wise.github.io/probmed/reference/pmedW_md.md).

## Usage

``` r
sinkhorn_div(X, Y, eps = 0.1, iters = 2000L)
```

## Arguments

- X, Y:

  Numeric matrices (rows = samples, columns = dimensions).

- eps:

  Numeric: regularization parameter (default 0.1).

- iters:

  Integer: maximum Sinkhorn iterations (default 400).

## Value

Non-negative scalar.
