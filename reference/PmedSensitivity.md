# P_med Sensitivity Result (M–Y Confounding)

S7 result for the sensitivity of an indirect-to-total effect ratio to an
additive bias in its numerator. It stores the original ratio, the ratio
on the bias grid, and the bias needed to make the indirect effect zero
or to reach a specified ratio threshold.

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

  Numeric: the ratio target used for `tipping_threshold` (`NA` if not
  requested).

- evalue:

  Numeric: optional E-value summary (length-2: estimate and CI bound)
  from the EValue package; `numeric(0)` if not computed.

- scale:

  Character: label of the ratio scale (e.g. `"Delta"`, `"incremental"`);
  free text supplied by the caller.

- call:

  Call: the originating call.

## See also

[`pmed_sensitivity()`](https://data-wise.github.io/probmed/reference/pmed_sensitivity.md)
