# Gauge-Calibrated Proportion Mediated Result

S7 class for the gauge-calibrated proportion mediated. Reports the
interventional proportion mediated `P_med = IIE / OE` alongside the
**gauge residual** `W = R / OE`, where `R = OE - IDE - IIE` is the
non-decomposability (treatment-by-mediator interaction) term. A `W`
significantly different from zero signals that the additive split of the
overall effect fails and the single-number `P_med` is unreliable.

## Usage

``` r
GaugePmedResult(
  p_med = integer(0),
  p_med_ci = integer(0),
  p_med_fieller = numeric(0),
  fieller_type = NA_character_,
  W = integer(0),
  W_ci = integer(0),
  W_se = integer(0),
  W_p = integer(0),
  OE = integer(0),
  IDE = integer(0),
  IIE = integer(0),
  R = integer(0),
  theta = integer(0),
  method = character(0),
  n = integer(0),
  ci_level = integer(0),
  se_method = "analytic",
  reps = 1L,
  call = NULL
)
```

## Arguments

- p_med:

  Numeric: interventional proportion mediated, `IIE / OE`.

- p_med_ci:

  Numeric length-2: confidence interval for `p_med`.

- p_med_fieller:

  Numeric: Fieller confidence set for the ratio `p_med` (may be
  unbounded; empty if not requested).

- fieller_type:

  Character: type of Fieller set returned (e.g. bounded, unbounded, or
  empty); `NA` if not computed.

- W:

  Numeric: gauge residual `R / OE`.

- W_ci:

  Numeric length-2: confidence interval for `W`.

- W_se:

  Numeric: standard error of `W`.

- W_p:

  Numeric: two-sided p-value for `H0: W = 0`.

- OE, IDE, IIE, R:

  Numeric: overall, interventional direct, interventional indirect
  effects and the remainder.

- theta:

  Numeric length-4: corner means `theta(a, a')`.

- method:

  Character: estimation method.

- n:

  Integer: sample size.

- ci_level:

  Numeric: confidence level.

- se_method:

  Character: `"analytic"` (influence-function se) or `"bootstrap"`
  (nonparametric, valid but mildly conservative).

- reps:

  Integer: number of repeated cross-fitting fold draws averaged for the
  point estimate.

- call:

  Call: original call.
