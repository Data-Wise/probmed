# Renormalization-Group Flow of the Proportion Mediated

S7 class for the **scale-indexed proportion mediated**
`g_l = NIE_l / TE_l` estimated on data coarse-grained to aggregation
scale `l` (block size `m`). Under the coarse-graining operator `T_m` –
which replaces each within-cluster block of `m` units by its mean –
`g_l` flows with `m`, and the flow is governed by the mediation
beta-function.

This is the classical **proportion mediated** (a ratio of effects),
*not* the package's `P_med` (a probability in `[0, 1]`); the two are
different quantities and are deliberately named apart. See
`vignette("probmed")` for `P_med`.

## Usage

``` r
RgFlowResult(
  g = integer(0),
  g_ci = integer(0),
  g_se = integer(0),
  a_path = integer(0),
  b_path = integer(0),
  TE = integer(0),
  grad_norm = integer(0),
  phi_cluster = integer(0),
  m = integer(0),
  n_units = integer(0),
  n_clusters = integer(0),
  ci_level = integer(0),
  call = NULL
)
```

## Arguments

- g:

  Numeric: estimated proportion mediated at this scale.

- g_ci:

  Numeric length-2: cluster-robust EIF Wald confidence interval.

- g_se:

  Numeric: cluster-robust standard error.

- a_path, b_path, TE:

  Numeric: the `a`, `b`, and total-effect components.

- grad_norm:

  Numeric: Euclidean norm of the gradient of `g` in the five moments;
  bounds the bias induced by a moment perturbation of a given size.

- phi_cluster:

  Numeric: per-cluster influence contributions, in cluster order.
  Retained so two scales estimated on the same clusters can be
  contrasted by
  [`rg_flow_contrast()`](https://data-wise.github.io/probmed/dev/reference/rg_flow_contrast.md).

- m:

  Numeric: aggregation scale (block size); `1` is the unit scale.

- n_units:

  Integer: coarse-grained units retained at this scale.

- n_clusters:

  Integer: clusters contributing at this scale.

- ci_level:

  Numeric: confidence level.

- call:

  Call: original call.

## Details

`g_l` is a smooth functional of five second moments of `(A, M, Y)`, so
the efficient influence function is available in closed form and the
plug-in estimator is already efficient in the linear-Gaussian case (no
cross-fitting required). The sampling unit is the **cluster**, not the
unit: influence contributions are aggregated to the cluster before their
variance is taken.
