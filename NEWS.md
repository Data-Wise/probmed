# probmed 0.3.0.9000 (development)

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
