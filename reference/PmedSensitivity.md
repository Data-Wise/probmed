# P_med Sensitivity Result (M–Y Confounding)

S7 class holding a one-row sensitivity analysis for a
proportion-mediated estimand under violation of assumption **A4** (no
exposure-induced mediator–outcome confounding). Every estimand in the
`P_med` family shares the ratio form `P_med = indirect / total`. A4
failure biases the **cross-world** corner mean that enters the
*indirect* (numerator) component, while the *total* (denominator) effect
is identified from observable corner means and is robust to A4. The
object therefore reports how `P_med` moves as the numerator bias `b`
ranges over a stated set, plus the **tipping bias** that drives the
indirect component (or `P_med`) through a threshold.

## Usage

``` r
PmedSensitivity(
  p_med = integer(0),
  indirect = integer(0),
  total = integer(0),
  bias_grid = integer(0),
  p_med_bias = integer(0),
  tipping_indirect = integer(0),
  tipping_threshold = NA_real_,
  threshold = NA_real_,
  evalue = numeric(0),
  scale = NA_character_,
  call = NULL
)
```

## Arguments

- p_med:

  Numeric: the point estimate `indirect / total` (b = 0).

- indirect:

  Numeric: the indirect (mediated) component, the numerator.

- total:

  Numeric: the total effect, the denominator.

- bias_grid:

  Numeric: the grid of numerator biases `b` evaluated.

- p_med_bias:

  Numeric: `P_med` at each `b`, i.e. `(indirect + b)/total`.

- tipping_indirect:

  Numeric: bias `b` at which `indirect + b = 0` (the indirect effect
  vanishes). `NA` if `total == 0`.

- tipping_threshold:

  Numeric: bias `b` at which `P_med` hits `threshold`; `NA` if not
  requested.

- threshold:

  Numeric: the `P_med` target used for `tipping_threshold` (`NA` if not
  requested).

- evalue:

  Numeric: optional E-value summary (length-2: estimate and CI bound)
  from the EValue package; `numeric(0)` if not computed.

- scale:

  Character: label of the `P_med` scale (e.g. `"Delta"`,
  `"incremental"`); free text supplied by the caller.

- call:

  Call: the originating call.

## See also

[`pmed_sensitivity()`](https://data-wise.github.io/probmed/reference/pmed_sensitivity.md)
