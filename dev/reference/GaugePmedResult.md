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
  W_ci_wald = numeric(0),
  weak_id = NA,
  weak_id_ratio = NA_real_,
  oe_snr = NA_real_,
  oe_regular = NA,
  weak_id_ratio_threshold = 3,
  oe_snr_threshold = 2,
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

- W_ci_wald:

  Numeric length-2: symmetric Wald interval for `W`, retained as the
  reference for the weak-identification diagnostic even under
  `se_method = "bootstrap"` (where `W_ci` holds the percentile
  interval).

- weak_id:

  Logical: `TRUE` when the percentile CI for `W` is at least
  `weak_id_ratio_threshold` (default 3) times wider than the symmetric
  Wald interval. `NA` unless `se_method = "bootstrap"` (both intervals
  are required), and `NA` (not `FALSE`) when the Wald interval has zero
  width. The ratio is a scale-only adaptation of the
  bootstrap-vs-asymptotic- discrepancy idea in Zhan (2026); it is not
  Zhan's statistic, which is a Kolmogorov-Smirnov distance for linear
  IV.

  **What it measures** (2026-08-22, corrected corner weights;
  `docs/specs/FINDINGS-2026-08-22-postfix-bootstrap-gates.md`): a
  bootstrap-vs-Wald width discrepancy, a symptom of bootstrap
  instability near `OE = 0` – **not a coverage diagnostic**. Against
  exact truth (4 cells x 200 datasets, n = 800) draws flagged at *any*
  threshold from 1.2 to 3 have the same Wald coverage as unflagged draws
  (0.94-0.95 vs 0.94-0.97). In the regular cells the ratio is 1.02
  (median), 1.22 (q95), 1.49 (q99); near the null it tracks
  `oe_regular = FALSE` (P(flag \| irregular) 0.82-0.88 at thresholds \>=
  1.5 vs 0.15-0.40 given regular) and at 3x adds 2 of the 98 near-null
  draws that pass A1. No
  [`warning()`](https://rdrr.io/r/base/warning.html) is issued on it;
  [`print()`](https://rdrr.io/r/base/print.html) reports the ratio when
  it exceeds the threshold. `oe_regular` is the gate to consult, with
  the Fieller set.

  **History.** The flag shipped with its 3x threshold "calibrated above
  the ~2x width gap of the anti-conservative Wald interval" and was
  validated on a 48,000-rep grid as a "narrowness proxy" whose coverage
  separation was explained by Wald width
  (`inst/sim/results/weakid_validation_cells.csv`). All of that measured
  the estimator *before* the corner-EIF weight fix (NEWS 0.3.0.9000):
  the ~30% se shortfall and the 2x gap were the bug. The threshold stays
  at 3 – with no signal to calibrate to, a new number would imply a
  calibration that does not exist – and the field is kept because a user
  running the bootstrap may want to see its instability.

- weak_id_ratio:

  Numeric: ratio of the percentile `W`-interval width to the Wald
  `W`-interval width; `NA` if not computed.

- oe_snr:

  Numeric: signal-to-noise of the denominator, `|OE| / se(OE)`. Small
  values indicate `OE` near 0, where `W = R/OE` is non-regular.

- oe_regular:

  Logical: `TRUE` when `oe_snr >= 2`; `FALSE` flags a near-singular `OE`
  – the denominator of the ratio `W = R/OE`, whose interval estimation
  is the classical Fieller (1954) problem. Computed unconditionally from
  the analytic se(OE) regardless of `se_method` (only the accompanying
  [`warning()`](https://rdrr.io/r/base/warning.html) is gated on
  `se_method = "bootstrap"`). `NA` only in the degenerate case
  `se(OE) == 0`.

  **What `FALSE` means – width, not under-coverage.** (Re-measured with
  the corrected corner weights, 2026-08-22, same reading: in the
  near-null cell flagged draws cover 0.99 at a median Wald width 8.7x
  \|W\|, while the near-null draws that *pass* the gate cover 0.90 at
  4.1x – `docs/specs/FINDINGS-2026-08-22-postfix-bootstrap-gates.md`.)
  In the 48,000-rep validation grid, draws flagged by this gate did
  **not** under-cover; they *over*-covered (+0.10, in 22/22 qualifying
  cells, stable across every design slice). The mechanism: as `OE`
  approaches 0 the ratio `W = R/OE` and its intervals blow up – flagged
  Wald intervals were a median **21x wider than the true \|W\|** (vs 3x
  unflagged). Such an interval contains the truth and nearly everything
  else; it is untrustworthy because it is **uninformative**, not because
  it misses. Of the two gates this is the stronger diagnostic – its
  effect survives controlling for interval width (+0.036, SE 0.003),
  unlike `weak_id`'s. Treat `oe_regular = FALSE` as "the point estimate
  and CI for `W` are numerically meaningless here", not as a coverage
  warning.

- weak_id_ratio_threshold, oe_snr_threshold:

  Numeric: the threshold values actually used for the
  `weak_id`/`oe_regular` gates on this result (see the identically-named
  arguments of
  [`ward_residual()`](https://data-wise.github.io/probmed/dev/reference/ward_residual.md)).
  Stored so [`print()`](https://rdrr.io/r/base/print.html) can report
  the threshold that was actually applied, not a hardcoded default.

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
