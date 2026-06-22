# Incremental Mediated Elasticity Result

S7 class for the **incremental mediated elasticity**
`P_med^delta(delta)` (Paper 2 of the P_med-for-modern-estimands
program). Rather than a single number, the estimand is a *curve*: the
share of the marginal (derivative-scale) total effect that flows through
the mediator when the treatment-assignment odds are tilted by a factor
`delta`.

Writing `Theta(delta_1, delta_2)` for the mean outcome under a
`delta_1`-tilt of the direct pathway and a `delta_2`-tilt of the
mediated pathway, the elasticities are `theta'_dir = d/d delta_1 Theta`,
`theta'_med = d/d delta_2 Theta`, and
`theta'_tot = d/d delta Theta(delta, delta)`. The multivariate chain
rule gives the **exact decomposition**
`theta'_tot = theta'_dir + theta'_med` (remainder `R(delta) = 0` for
every `delta`), and `P_med^delta(delta) = theta'_med / theta'_tot`.

Under a linear-Gaussian model with no `A * M` interaction the curve is
flat and equals the classical proportion mediated; it bends with `delta`
only when an interaction is present.

## Usage

``` r
IncrPmedResult(
  curve = data.frame(),
  deltas = integer(0),
  method = character(0),
  n = integer(0),
  ci_level = integer(0),
  call = NULL
)
```

## Arguments

- curve:

  Data frame: one row per `delta`, columns `delta`, `dir`, `med`, `tot`,
  `Pmed`, `se`, `lo`, `hi`.

- deltas:

  Numeric: the tilt factors at which the curve was evaluated.

- method:

  Character: estimation method.

- n:

  Integer: sample size.

- ci_level:

  Numeric: confidence level.

- call:

  Call: original call.
