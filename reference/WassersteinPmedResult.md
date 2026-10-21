# Wasserstein / Transport-Scale Proportion Mediated Result

S7 class for the **Wasserstein proportion mediated**
`P_med^W = NIE^W / (NIE^W + NDE^W)`, where `NIE^W = W_2(nu_11, nu_10)`
and `NDE^W = W_2(nu_10, nu_00)` are the 2-Wasserstein distances between
the corner laws `nu(a, a') = Law(Y(a, M(a')))`. Unlike the additive
split `TE = NDE + NIE`, the Wasserstein decomposition handles
treatment-by-mediator interaction via the synergy term
`G^W = NIE^W + NDE^W - TE^W >= 0` (triangle gap).

## Usage

``` r
WassersteinPmedResult(
  pmedW = integer(0),
  pmedW_ci = integer(0),
  pmedW_se = integer(0),
  NIE_W = integer(0),
  NDE_W = integer(0),
  TE_W = integer(0),
  synergy = integer(0),
  boundary = FALSE,
  method = character(0),
  n = integer(0),
  ci_level = integer(0),
  call = NULL
)
```

## Arguments

- pmedW:

  Numeric: transport-scale proportion mediated.

- pmedW_ci:

  Numeric length-2: confidence interval (Wald or bootstrap).

- pmedW_se:

  Numeric: standard error.

- NIE_W, NDE_W, TE_W:

  Numeric: Wasserstein NIE, NDE, and total effect.

- synergy:

  Numeric: triangle gap `NIE^W + NDE^W - TE^W >= 0`.

- boundary:

  Logical: `TRUE` when either Wasserstein component is near zero.

- method:

  Character: estimation method (`"eif"`, `"bootstrap"`, or `"dr"`).

- n:

  Integer: sample size.

- ci_level:

  Numeric: confidence level.

- call:

  Call: original call.
