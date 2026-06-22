# Sobol / Variance-Scale Proportion Mediated Result

S7 class for the **Sobol (functional-ANOVA) proportion mediated**
`P_med^{sigma2} = V_med / V_T` (Paper 3 of the
P_med-for-modern-estimands program). On the variance scale the
proportion mediated is the first-order Sobol sensitivity index of the
outcome with respect to the mediated pathway: the share of the total
interventional outcome variance `V_T` explained by the mediated
component `V_med`. Unlike the additive split `TE = NDE + NIE`, the
variance decomposition `V_T = V_dir + V_med + V_int` always holds (the
interaction variance `V_int` is a separate non-negative term), so
`P_med^{sigma2}` is well defined even under treatment-by-mediator
interaction.

## Usage

``` r
SobolPmedResult(
  p_med = integer(0),
  se = integer(0),
  ci = integer(0),
  ci_wald = integer(0),
  ci_A = NA_real_,
  ci_B1 = NA_real_,
  ci_B2 = NA_real_,
  procedure = "B",
  boundary = logical(0),
  Pmed_upper = integer(0),
  S1_med = integer(0),
  ST_med = integer(0),
  Vd = integer(0),
  Vm = integer(0),
  Vdm = integer(0),
  VT = integer(0),
  Dm = integer(0),
  se_Dm = integer(0),
  vmed_split_p = NA_real_,
  vmed_split_reject = NA_integer_,
  theta = integer(0),
  method = character(0),
  se_method = "analytic",
  reps = 1L,
  n = integer(0),
  ci_level = integer(0),
  call = NULL
)
```

## Arguments

- p_med:

  Numeric: Sobol proportion mediated `V_med / V_T`.

- se:

  Numeric: standard error of `p_med` (delta-method, ratio identity).

- ci:

  Numeric length-2: **reported** confidence interval, selected by
  `procedure` – Procedure B (default, uniformly valid) or Procedure A
  (gated).

- ci_wald:

  Numeric length-2: ungated symmetric Wald interval (always the Wald
  interval, even at the boundary).

- ci_A:

  Numeric length-2: Procedure-A (gated, pre-test) interval – the Wald
  interval off-boundary, or `[0, Pmed_upper]` at the boundary.
  Pointwise- but **not** uniformly-valid (see Coverage caveat).

- ci_B1:

  Numeric length-2: Procedure-B interval from the **same-sample**
  `Delta_m` Wald CI mapped through `delta -> (c_m/V_T) delta^2`
  (efficient).

- ci_B2:

  Numeric length-2: **experimental** Procedure-B variant from a single
  decorrelated split (`Delta_m` point on one half, se on the independent
  half) mapped through the same square. Intended to break `Delta_m`/se
  coupling, but the half-sample centre makes it high-variance and it did
  **not** improve near-null coverage in simulation – a
  K-fold-decorrelated se keeping the full-sample point is the open fix.
  Not recommended; use `ci` (B1).

- procedure:

  Character: `"B"` (default, uniform sample-split CI) or `"A"` (legacy
  gated). Controls which interval is returned in `ci`.

- boundary:

  Logical: `TRUE` when the split test does not reject `H0: V_med = 0`
  (the symmetric Wald CI is non-regular there).

- Pmed_upper:

  Numeric: one-sided Procedure-A upper bound for `p_med`.

- S1_med:

  Numeric: first-order Sobol index `V_med / V_T` (equals `p_med`).

- ST_med:

  Numeric: total Sobol index `(V_med + V_int) / V_T`.

- Vd, Vm, Vdm, VT:

  Numeric: direct, mediated, interaction, and total variance components.

- Dm:

  Numeric: mediated antecedent contrast `Delta_m`.

- se_Dm:

  Numeric: standard error of `Delta_m`.

- vmed_split_p:

  Numeric: p-value of the sample-split boundary test of `H0: V_med = 0`.

- vmed_split_reject:

  Integer: `1L` if the boundary test rejects, `0L` otherwise (`NA` when
  the test was not run).

- theta:

  Numeric length-4: corner means `theta(a, a')`.

- method:

  Character: estimation method.

- se_method:

  Character: `"analytic"` (default, influence-function se) or
  `"bootstrap"` (nonparametric resample-and-refit se; valid but
  conservative near the boundary). See the boundary coverage section.

- reps:

  Integer: number of repeated cross-fitting fold draws averaged for the
  point estimate (default `1`); `reps > 1` removes the near-boundary
  fold-split Monte-Carlo variance (mean DML aggregation).

- n:

  Integer: sample size.

- ci_level:

  Numeric: confidence level.

- call:

  Call: original call.

## Details

Writing the four corner means as `theta(a, a') = E[Y(a, M(a'))]`, the
contrasts are the direct antecedent `Delta_d`, the mediated antecedent
`Delta_m`, and the interaction remainder
`R = theta_11 - theta_10 - theta_01 + theta_00`. With Bernoulli tilt
variances `c_d = p_d(1 - p_d)` and `c_m = p_m(1 - p_m)`,
`V_dir = c_d Delta_d^2`, `V_med = c_m Delta_m^2`, `V_int = c_d c_m R^2`,
and `V_T = V_dir + V_med + V_int`. The estimator is the cross-fitted
one-step estimator built on the shared triply-robust corner
pseudo-outcomes (the same engine as
[`ward_residual()`](https://data-wise.github.io/probmed/reference/ward_residual.md)
and
[`incr_pmed()`](https://data-wise.github.io/probmed/reference/incr_pmed.md));
the delta method propagates the corner influence functions through the
quadratic `V`-components and the ratio identity to a standard error for
`P_med^{sigma2}`.

**Boundary remedy.** Because `V_med = c_m Delta_m^2` is a squared
quantity, its influence function
`phi_{V_med} = 2 c_m Delta_m phi_{Delta_m}` degenerates as
`Delta_m -> 0`, deflating the Wald standard error and collapsing
coverage at the null `V_med = 0`. The null is equivalent to
`Delta_m = 0`, which has a regular root-`n` influence function.
Following Williamson et al. (2021, *Biometrics*), the boundary is tested
with a sample-split statistic that estimates `Delta_m` on one half and
its standard error on the independent other half (so `z ~ N(0, 1)` even
at the boundary).

**Two interval procedures.** The default reported `ci` is **Procedure
B**: the image of the regular two-sided `Delta_m` Wald CI under the map
`delta -> (c_m/V_T) delta^2` (`ci_B1`/`ci_B2`). Because it is the
continuous image of a uniformly-valid CI for the *regular* parameter
`Delta_m`, it is uniformly valid in `Delta_m` – there is no pre-test,
hence no boundary pathology (subject to the near-null se caveat below).
The reported `ci` is `ci_B1`, the same-sample image (centered at the
full `Delta_m`, so it always brackets `p_med`); `ci_B2` is an
*experimental* single-split variant that did not improve coverage (see
its note). **Procedure A** (`procedure = "A"`, `ci_A`) is the legacy
*gated* rule: Wald on rejection, one-sided `[0, Pmed_upper]` otherwise.

**Coverage near the boundary – mechanism and remedies (A-15).** Two
distinct issues; the second is now decomposed (reproducible fixed-seed
variance decomposition) into three findings, each with a concrete
remedy.

*(1) Procedure A is not uniformly valid.* As a pre-test (gating) rule A
is only pointwise valid: near the boundary it routes downward-selected
`Delta_m_hat` to the contracting one-sided bound, a Leeb-Potscher
post-selection effect. Prefer B, which does not gate.

*(2) The near-null `Delta_m` standard error.* Three pinned facts:

- **Fold-split Monte-Carlo variance dominates.** A single cross-fitting
  fold partition contributes ~80% of `Var(Delta_m_hat)` in the near-null
  regime (~45% at ordinary effect sizes), persistent across `n = 2000`
  to `8000`. This is an algorithmic nuisance, not sampling information:
  it inflates the estimator's variance and makes the single-split point
  estimate seed-dependent. *Remedy:* **`reps > 1`** (repeated
  cross-fitting – the corner influence matrix is averaged over `reps`
  independent fold draws), which removes it and yields a reproducible,
  lower-variance point estimate (mean DML aggregation).

- **`Delta_m_hat` is approximately normal at the null.** Recomputing
  coverage with the oracle Monte-Carlo SD gives ~0.95, so the Wald
  *shape* is correct for the regular contrast `Delta_m`: the
  under-coverage is a wrong interval *width*, not a wrong *shape* (there
  is no boundary non-normality for `Delta_m` itself).

- **A residual analytic-se bias remains.** Even after fold noise is
  removed, the analytic influence-function se for `Delta_m` is ~0.8x its
  true sampling SD near the null (structural; persistent in `n`); it is
  well calibrated (ratio ~1) at ordinary effect sizes. So the *default*
  analytic near-null interval – and its Procedure-B image – covers only
  ~0.85. *Remedy:* **`se_method = "bootstrap"`** (resample rows, refit
  the cross-fit estimator), which recovers a valid but mildly
  **conservative** se (~1.25x the true SD, coverage ~0.97 near null) and
  is calibrated off-boundary.

**Practical guidance.** Off the boundary the analytic default is
calibrated and fast. Near the boundary (the split test does not reject),
pass `reps > 1` for a reproducible point and `se_method = "bootstrap"`
for a valid interval (conservative there – it trades width for
guaranteed coverage). `V_T` is well calibrated throughout. The full
coverage grid across the `Delta_m` transition is validated by a separate
large simulation. Use `procedure = "A"` only to reproduce the legacy
gated behaviour.

## References

Williamson, B. D., Gilbert, P. B., Carone, M., & Simon, N. (2021).
Nonparametric variable importance assessment using machine learning
techniques. *Biometrics*, 77(1), 9–22.
