# probmed 0.3.0.9000 (development)

## `weak_id` no longer warns; bootstrap arm and gates re-measured with correct weights

* Re-measured after the weight fix below (`inst/sim/boot_gates_check.R`,
  `inst/sim/boot_gates_calibrate.R`, `inst/sim/results/boot_gates_postfix.csv`;
  `docs/specs/FINDINGS-2026-08-22-postfix-bootstrap-gates.md`; 4 cells × 200
  datasets × B = 200 against exact truth):
  - **Bootstrap (percentile) arm:** coverage 0.905 / 0.930 / 0.920 in the
    continuous / binary / intermediate cells (Wald 0.935 / 0.970 / 0.910) at the
    same width as Wald; bootstrap SD unstable for binary `Y` (CV 2.2). No longer
    "~1.00"; never better than the default. Kept as an option.
  - **A2 `weak_id` cannot be calibrated to coverage:** at every threshold from
    1.2 to 3, flagged draws cover as well as unflagged (0.94–0.95 vs 0.94–0.97).
    The width ratio tracks bootstrap instability, which `oe_regular` already
    flags; beyond A1 it catches 2 of 98 near-null draws. `ward_residual()`
    **no longer issues a warning** on `weak_id`; the field, ratio and threshold
    (still 3) are kept and `print()` reports the ratio as a "bootstrap/Wald
    width discrepancy". Docs rewritten (`weak_id`, `weak_id_ratio_threshold`,
    `se_method`, vignette).
  - **A1 `oe_regular` stands:** flagged near-null draws over-cover (0.99 at a
    median Wald width 8.7× |W|); near-null draws that pass it cover 0.90 — the
    local-to-zero regime no sample gate separates; use the Fieller set.

## Bug fix: corner-EIF inverse-probability weights were one constant per fold

* `.corner_fit()` — the cross-fit corner influence engine behind
  `ward_residual()`, `incr_pmed()` and `sobol_pmed()` — built its weights as
  `pa <- function(z) ifelse(z == 1, p1, 1 - p1)`. `ifelse()` returns a result
  shaped like its **test**, and `z` is a scalar, so `pa(a)` was `p1[1]`: the
  propensity of the **first test row** of the fold, applied to every row. The
  mediator-density proxy `qa()` had the same defect. Shipped in `10b41aa`
  (PR #8), copied into `pmedW_dr()` (`R/wasserstein-pmed.R`) and into the
  HPC standalone `inst/sim/gauge_coverage_standalone.R`; all three fixed.
  Point estimates were not visibly biased (the EIF is triply robust and the
  outcome / projection models were correctly specified) but the influence
  function was not the efficient one and its variance depended on which row
  happened to land first in each fold.

* **What the fix changes, measured** (`inst/sim/se_shipped_check.R`, exact
  truth, 250 datasets per cell, n = 800; numbers are median se/empSD, CV of
  the se across datasets, Wald coverage):

  | cell | `reps = 1` | `reps = 10` | oracle |
  |---|---|---|---|
  | continuous, strong ID | 0.89 / 0.23 / **0.936** | 0.91 / 0.22 / **0.936** | 0.87–0.90 / 0.18 / 0.940 |
  | binary, strong ID | 0.99 / 0.33 / **0.972** | 0.97 / 0.32 / **0.968** | 0.95 / 0.31 / 0.972 |
  | intermediate (s = 0.5, τ = 0.4) | 0.95 / 0.29 / 0.912 | — | 0.90 / 0.25 / 0.964 |
  | near-null (s = 0.2) | non-regular (CV 3.4) | — | non-regular (CV 8.8) |

  The single-partition estimator's empSD is now 0.063 vs oracle 0.056
  (continuous) and 0.121 vs 0.122 (binary) — before the fix it was 0.121 and
  0.175. The default analytic Wald interval is calibrated in the regular
  cells; `reps > 1` buys little; near the null the ratio is non-regular for
  every estimator, as before.

* **Superseded.** Every coverage figure the package reported before this
  entry was produced with the bug: the 2,000-rep grid
  (`gauge_boot_coverage_nsim2000.csv`: Wald ~0.85–0.90, percentile ~1.00), the
  48,000-rep `weak_id`/`oe_regular` validation grid, and PR #33's variance
  decomposition (next section). The "fold-split noise" those diagnosed was the
  bug — a different first row per partition is a different constant weight per
  partition. The bootstrap arm and both gates have **not** been re-measured
  with the fix; their documentation says so. A full-sample-nuisance se for
  `reps > 1` was built on this branch as the "fix" for the apparent se
  instability, measured, and dropped once the real cause was found
  (`.corner_fit_full()` stays as a deterministic test probe;
  `inst/sim/se_candidates{,2,3}.R`, `se_bootstrap_reps.R` record the search).

* Tests: row-order invariance of the full-sample corner means (fails on the
  old weights), pinned post-fix values, and the `reps > 1` se construction
  rebuilt by hand.

* **Second fix, exposed by the first: `incr_pmed()`'s g-score term.** Term II
  of the influence function multiplied `A - g(C)` by the **per-row**
  `a_med`/`a_dir`, which carry the corner EIF's inverse-probability noise; the
  orthogonality correction needs their conditional means given `C`. With the
  old constant-per-fold weights this looked calibrated (se/empSD
  1.07 / 1.11 / 1.14 at δ = 0.5 / 1 / 2, inside the test's band); with correct
  weights it inflated to 0.97 / 1.15 / **1.40**. Now projected linearly on the
  covariates (the same device as `eta` in `.corner_phi()`): 0.94 / 0.95 / 0.97.
  The weight fix alone also made `incr_pmed()` ~1.5x more efficient (empSD of
  `P_med^delta` 0.083 → 0.055 at n = 1000). Point estimates unchanged by the
  projection (Term II is mean-zero).

## Simulation findings (PR #33) — superseded by the weight fix above

* **Where `ward_residual()`'s variance comes from, and why its analytic se
  under-covers** — `inst/sim/phi_decomposition.R`, results in
  `inst/sim/results/phi_decomp_summary.csv`, full account in
  `docs/specs/FINDINGS-2026-08-22-phi-decomposition.md`. The coverage grid's
  `se/empSD = 0.65-0.83` was a **skewness artifact of an unstable se**, not a
  biased formula: with the DGP's true nuisances plugged into the same corner EIF
  the se is exact and coverage nominal. In the grid's regular cells the
  single-partition estimator's excess variance is **fold-split noise**
  (`.corner_fit()` redraws the partition every call); `reps > 1`, already
  shipped, brings the point estimate to oracle efficiency for continuous `Y`
  (n=800: empSD 0.121 → 0.057 vs oracle 0.056). Coverage stays ~0.88 only
  because the se estimator has CV 0.5–0.8 across datasets (oracle 0.17) — a
  better se for the `reps > 1` estimator was the open item. [All of this
  measured the pre-fix weights; see the bug-fix entry above.] Near the null even
  the oracle se explodes: the ratio is non-regular there, which is exactly what
  `oe_regular` flags. Binary `Y` with a small `OE` is a transition regime where
  `reps = 10` halves the variance but does not reach the oracle.

* Docs corrected accordingly (`se_method`, `weak_id`, `oe_regular`; the
  gauge-residual vignette). Two citation fixes: the `weak_id` gate is a
  scale-only **adaptation** of Zhan (2026)'s bootstrap-vs-asymptotic idea, not
  his statistic (a Kolmogorov–Smirnov test for linear IV); and the
  bootstrap-consistency result of Lin and Han (2026) holds the nuisances fixed
  and does **not** cover the package's refit-per-resample bootstrap, which the
  docs had claimed. The user-facing warnings no longer cite either.

* `AGENTS.md` and the `gauge_sim_out` scratch directory are now
  `.Rbuildignore`d (they raised a top-level-files NOTE in `R CMD check`).

## Documentation

* `ward_residual()` / `GaugePmedResult`: document what the `weak_id` flag does
  and does **not** tell you. Two things a user needs and the docs did not say:
  (1) **an absent flag does not buy nominal coverage** — independently of the
  flag, the Wald interval for `W` is anti-conservative everywhere (~0.88 across
  the 16,000-rep coverage grid, never nominal in any cell), so
  `weak_id = FALSE` at best means "not the worst tail of an already sub-nominal
  default". `se_method = "bootstrap"` is the safer arm, but know what it costs:
  the percentile interval covered **1.00** in all 8 cells — uncalibrated in the
  safe direction, i.e. wide enough to be uninformative near the null. **Neither
  arm is nominal**; pick the error you can live with. (2) The flag trades
  sensitivity for specificity — but read the magnitudes as indications, not
  measurements. A 120-draw pilot (`inst/sim/pilot/weak_id_pilot.R`), conditioning
  per draw, gives false alarms `flag | oe_snr >= 3` of **0/58** (exact 95% CI
  `[0.00, 0.06]` — bounded near 6%, *not* shown to be zero) and detection
  `flag | oe_snr <= 1.2` of **12/28 = 0.43** (CI `[0.24, 0.63]`), with a further
  34 draws in neither regime. The *mechanism* is firmer than the magnitudes —
  `weak_id_ratio`'s draw-to-draw SD approaches its mean in the weak regime,
  which necessarily costs detection.

* The **#11 validation grid ran** (48,000 reps, 24 cells,
  `inst/sim/hopper/run_weakid_validation.R`; results in `inst/sim/results/`,
  adversarially re-verified from the raw reps) and settled the two gates'
  division of labour — differently from what the docs had anticipated:
  - `weak_id` (percentile/Wald width ratio ≥ 3): flagged draws do have worse
    Wald coverage in **all 24 cells**, but that separation is **almost entirely
    explained by the Wald interval's own width** (the flag's effect collapses
    to ~0 under a flexible width control; the bootstrap numerator adds nothing
    detectable). The docs now present it as a **self-contained narrowness
    proxy** — useful because a single fit offers no reference for "narrow" —
    not as evidence the bootstrap-vs-Wald comparison detects a distinct
    pathology. The threshold 3 is a **convention**: 1.75–4 were not
    distinguishable, and the grid spans one simulation design.
  - `oe_regular` (`oe_snr < 2`): its flagged draws ***over*-cover** (+0.10, in
    22/22 qualifying cells) — as `OE → 0` the intervals blow up (median 21×
    wider than the true `|W|`) and cover everything. The docs now state it
    flags **uninformatively wide** intervals, not under-coverage — and note it
    is the *stronger* diagnostic, retaining an effect beyond width alone
    (unlike `weak_id`).

## New features

* `rg_flow_contrast()` and the `RgFlowResult` class (PR #30) add a
  **scale-indexed proportion mediated** — the classical NIE/TE ratio tracked
  across coarse-graining scales (`rg_coarse_grain()`), with **cluster-robust
  EIF inference** for clustered data. Deliberately *not* a P_med variant:
  probmed's first multilevel estimator. Takes raw clustered data plus cluster
  labels; whether `medfit::MediationData` should grow a cluster slot is an open
  data-contract question.

* New data set `multilevel_designs` (the package's first shipped data): design
  metadata for three clustered study designs -- ECLS-K:1998-99, ECLS-K:2011
  (kindergarten class of 2010-11), and the public `mediation::student`
  teaching data -- so users planning a
  clustered mediation analysis can start from documented cluster counts rather
  than round numbers. Records sampled and participating cluster counts
  separately, since the participating count (the one that describes the
  analyzable data) is well below the sampled count in both ECLS-K cohorts.
  Every value is transcribed from the public source cited in the row; the
  reproducible build script lives in `data-raw/multilevel_designs.R`.

* `wasserstein_pmed()` and the `WassersteinPmedResult` class add the
  **Wasserstein / transport-scale proportion mediated** (`P_med^W`), an
  optimal-transport effect size decomposing the total transport cost into
  natural indirect (`NIE^W`), natural direct (`NDE^W`), and synergy (`G^W`)
  components. Supports the 1-D quantile-coupling path, an entropic-OT
  (Sinkhorn) path for multivariate mediators (`d > 1`), and a doubly-robust
  corner-LAW estimator (`method = "dr"`) with an automatic bootstrap-CI
  fallback near the `NIE^W = 0` boundary (#16, #19).

* `ward_residual()` gains two identification diagnostics for the gauge residual
  `W = R/OE`, reported as new `GaugePmedResult` fields with matching `warning()`s:
  a **weak-identification flag** (`weak_id` / `weak_id_ratio`) that fires when the
  percentile CI for `W` is >= 3x wider than the symmetric Wald interval — a
  Wald-vs-percentile divergence diagnostic (Zhan 2026) — and a **regularity
  guard** (`oe_regular` / `oe_snr`) flagging a near-singular denominator
  (`|OE|/se(OE) < 2`), where the bootstrap CI for `W` is not valid
  (Lin et al. 2026). The Wald interval is retained as `W_ci_wald` (#11).

## Bug fixes

* `PmedResult` validator: the `ci_lower <= ci_upper` check was silently
  inert -- an S7 validator returns only its last expression's value, so with
  the checks written as separate `if`-blocks the bounds-order message was
  computed and discarded, and only the final `ci_level` check could error.
  The validator now accumulates all messages and returns them together, so a
  reversed confidence interval is correctly rejected. Found while closing a
  test-coverage gap on the core `PmedResult` print/summary/plot methods.

* `ward_residual()`: `print()` for the weak-ID / near-singular-OE diagnostics
  (added above) now reports the threshold actually applied to the result
  (`weak_id_ratio_threshold` / `oe_snr_threshold`, retained as new
  `GaugePmedResult` fields) instead of the hardcoded defaults ("3x" / "< 2"),
  which misreported the gate whenever a non-default threshold was passed to
  `ward_residual()`. `weak_id`/`oe_regular` also now stay `NA` (rather than a
  false-confident `FALSE`) in the degenerate case of a zero-width Wald interval
  or zero se(OE) (#11).

* `incr_pmed()`: the g-score term of the efficient influence function now
  carries the tilt-derivative (`q'`) weight, matching the point estimate.
  Previously it used the bare corner contrasts, inflating standard errors
  with `delta` and causing confidence intervals to over-cover (se_ratio up
  to ~4.6). CIs are now calibrated across `delta` (#20).

# probmed 0.3.0

## New features

* `ward_residual()` and the `GaugePmedResult` class add the **gauge-calibrated
  proportion mediated** for interventional/stochastic effects. Alongside the
  interventional proportion mediated `P_med = IIE/OE`, it reports the **gauge
  residual** `W = R/OE` with `R = OE - IDE - IIE`, the treatment-by-mediator
  interaction (non-decomposability) term. A nonzero `W` flags that the additive
  split fails and the single-number `P_med` is unreliable. Cross-fitted one-step
  estimator built on the triply-robust EIFs of the corner means
  `theta(a,a') = E[Y(a, M(a'))]`. Inference for the skewed ratios `W` and `P_med`:
  analytic sqrt(n) Wald intervals by default, or tail-aware **percentile**
  intervals via `se_method = "bootstrap"`; repeated cross-fitting (`reps`) to
  remove fold-split variance; and a **Fieller** confidence set that is honest
  about being unbounded when the total effect is not significant. Any two-level
  exposure coding is supported through `a0` (reference) / `a1` (comparison) — a
  factor, `{1,2}`, `{-1,1}` — with `>2` levels rejected rather than silently
  subsetted. General covariates; binary or continuous outcome. (feature/gauge-pmed
  + feature/gauge-bootstrap-se; companion manuscript: pmed-modern/01-gauge-pmed.)

* `incr_pmed()` and the `IncrPmedResult` class add the **incremental mediated
  elasticity** `P_med^delta(delta)` — the derivative-scale proportion mediated as
  a function of the treatment-tilt factor `delta`. Unlike a single number it is a
  *curve*; by the multivariate chain rule its direct and mediated elasticities sum
  to the total exactly (remainder zero for every `delta`), and the curve is flat at
  the classical `P_med` when there is no treatment-by-mediator interaction. Same
  cross-fitted one-step machinery as `ward_residual()`; ratio-identity SEs use the
  full efficient influence function including the g-score (propensity) term
  (Kennedy 2019, Cor. 2 term II) with the propensity cross-fit, making the
  estimator Neyman-orthogonal in the propensity score (efficient under
  all-nuisance convergence at `n^{-1/4}`). (feature/incremental-pmed; companion
  manuscript: pmed-modern/02-incremental-pmed.)

* `incr_sensitivity()` adds per-`delta` M–Y unmeasured-confounding sensitivity for
  an `IncrPmedResult`. For each tilt factor it reports the additive numerator bias
  that zeroes the incremental mediated share (`tipping = -med`) and, optionally,
  the bias driving `P_med^delta` to a user `threshold` (`threshold*tot - med`).
  A thin wrapper over `pmed_sensitivity()`; the additive-offset model is exact for
  the mean-based Paper 2 share. (feature/gauge-bootstrap-se; companion manuscript:
  pmed-modern/02-incremental-pmed, issue #4.)

* `sobol_pmed()`, `sobol_from_theta()` and the `SobolPmedResult` class add the
  **Sobol / functional-ANOVA variance share** `P_med^{sigma^2} = V_med / V_T` — the
  fraction of the intervention-induced outcome variance carried by the mediator
  pathway, with `V_med = c_m * Delta_m^2`. Same cross-fitted one-step corner-EIF
  machinery as `ward_residual()`. At the no-mediation boundary the variance share is
  non-regular (`V_med = c_m Delta_m^2` is degenerate at `Delta_m = 0`), so inference
  reduces to the regular contrast `Delta_m`. Two interval procedures: the default
  **Procedure B** (`procedure = "B"`, the image of the regular `Delta_m` Wald CI under
  the squared map — no pre-test) and the legacy gated **Procedure A**
  (`procedure = "A"`). Near-boundary inference adds two options (A-15): **`reps`**
  (repeated cross-fitting — averages the corner influence matrix over `reps` fold
  draws, removing the ~80% fold-split Monte-Carlo variance that dominates
  `Var(Delta_m_hat)` near the null and yielding a reproducible point estimate), and
  **`se_method = "bootstrap"`** (nonparametric resample-and-refit se — valid, mildly
  conservative near the non-regular boundary where the analytic influence-function se
  is ~0.8x anti-conservative). Defaults (`reps = 1`, `se_method = "analytic"`) are
  unchanged. `Delta_m_hat` is approximately normal at the null (oracle-SD coverage
  ~0.95), so the Wald shape is correct — the near-null issue is interval width, not
  shape; see `?SobolPmedResult`. (feature/sobol-pmed; companion manuscript:
  pmed-modern/03-sobol-pmed.)

# probmed 0.2.0 (2026-06-11)

## Features

* `pmed()` now accepts a `medfit::ParallelMediationData` object and computes the
  **joint** P_med for k parallel mediators — the probability the outcome is
  higher with all mediators at their treated levels than at their control
  levels: `Phi(delta * sum(a*b) / sqrt(2*sum(b^2*Vm) + 2*Vy))`, recovering the
  single-mediator formula at k = 1. All four methods are supported (`plugin`,
  `parametric_bootstrap`, `nonparametric_bootstrap`, `mbco`), with the total
  indirect effect `sum(a_j b_j)`. Gaussian outcome and mediators only. See
  `vignette("parallel-mediation")`.

* `pmed(..., method = "mbco")` adds a deterministic Model-Based Constrained
  Optimization interval (Tofighi & Kelley, 2020): a likelihood-ratio interval
  for both P_med and the indirect effect `a * b`, obtained by inverting the
  constrained-likelihood test rather than by resampling. Gaussian outcome and
  mediator (with covariates) and any contrast `x_ref != x_value`; seed-free and
  grid-resolution-independent. Non-Gaussian models still use the bootstrap
  methods.

## Fixes

* `print(PmedResult)` interpretation line now shows the mediation estimand
  `P(Y(1, M(1)) > Y(1, M(0)))` (manuscript Definition 1) instead of the stale
  direct-effect notation `P(Y_{X*, M_X} > Y_{X, M_X})`. The computed value was
  already correct; only the displayed notation was wrong. A plain-language gloss
  was added beneath it. README cached output updated to match.

# probmed 0.1.0 (2026-06-06)

First GitHub release (non-CRAN).

## Features

* `pmed()` computes P_med — a scale-free probabilistic effect size for causal
  mediation — from a formula or a `medfit::MediationData` object, with plugin,
  parametric-bootstrap, and nonparametric-bootstrap methods, alongside the
  indirect effect (`a * b`).

## Fixes

* `pmed()` now computes the **mediation** estimand
  `P(Y(x, M(x)) > Y(x, M(x*))) + 0.5 P(=)` (manuscript Definition 1): treatment
  held fixed, mediator varied between its treated and control levels with the
  tie term. Previously it computed a direct-effect contrast that depended on the
  direct effect `c'` and could land on the wrong side of 0.5. Verified against
  the closed form `Phi(ab / sqrt(2 (b^2 sigma_M^2 + sigma_Y^2)))` and the
  manuscript `memory_exp` example (`P_med = 0.68`).
* Binary/non-Gaussian outcomes now draw Bernoulli responses through the link
  (previously degenerate, returning 0/1), using the new medfit family slot.
* Parametric bootstrap indexes coefficients **by name** (was selecting
  intercepts as `a`/`b`); nonparametric bootstrap refits on the **correct
  family** (was always Gaussian).

## Ecosystem

* Builds on **medfit (>= 0.3.0)** for model extraction and the family/link slot.
  Part of the mediationverse ecosystem.

---

# probmed 0.0.0.9000

## Major Changes

* **S7 Architecture**: The package has been refactored to use the **S7** object-oriented system for robust class definitions and method dispatch.
* **Website Redesign**: The documentation website now uses the `litera` theme with a clean, academic design matching the `rmediation` package style.
* **Quarto Integration**: The package now uses Quarto (`.qmd`) for the README and vignettes, providing modern publishing capabilities.

## New Features

* **$P_{med}$ Calculation**: Implemented `pmed()` function to compute the probabilistic effect size $P_{med}$ for mediation analysis.
* **GLM Support**: Added support for Generalized Linear Models (e.g., logistic regression) for both mediator and outcome.
* **Bootstrap Inference**: Added parametric and nonparametric bootstrap methods for confidence intervals.
* **lavaan Integration**: Added `extract_mediation()` support for SEM models fitted with the `lavaan` package, including FIML and robust estimators.
* **mediation Integration**: Added support for extracting mediation structures directly from `mediation::mediate()` objects.
* **Indirect Effect Reporting**: `pmed()` now reports the Indirect Effect (product of coefficients) alongside $P_{med}$, including bootstrap confidence intervals.

## Documentation

* **Expanded README**: Added detailed explanation of $P_{med}$, features, and examples.
* **New Vignette**: Added "Introduction to probmed" vignette demonstrating linear and binary outcome examples.
* **Integration Vignettes**: Added dedicated vignettes for `lavaan` and `mediation` package integrations.
* **Comparison Vignette**: Added "Comparing probmed Workflows" to guide users on choosing the best integration method.
