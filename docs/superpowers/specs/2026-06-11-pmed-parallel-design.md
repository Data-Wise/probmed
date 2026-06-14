# Spec: Parallel (joint) P_med — `method` support for `ParallelMediationData`

**Date:** 2026-06-11
**Status:** approved (brainstorm) → ready for plan
**Package:** probmed
**Depends on:** medfit (`ParallelMediationData`, Gaussian extraction)

---

## 1. Summary

Extend P_med from a single mediator to **parallel multiple mediators**
(X → M₁, X → M₂, …, X → M_k, with all Mⱼ entering one outcome equation
Y ~ X + M₁ + … + M_k). This spec covers the **joint** estimand only: P_med for
moving *all* mediators together from their control to their treated levels.

The entry point is a new S7 method `pmed(<medfit::ParallelMediationData>, …)`
returning the existing `PmedResult`. All four inference methods are supported:
`plugin`, `parametric_bootstrap`, `nonparametric_bootstrap`, `mbco`.

**Scope decisions (locked in brainstorm):**

- **Parallel only** (serial/sequential is a separate later spec).
- **Joint estimand only** (path-specific per-mediator P_med deferred).
- **All four methods.**
- **Closed-form kernel (Approach A)** — exact Gaussian formula, no simulation.
- **Object-entry only** — `pmed(ParallelMediationData)`; a formula convenience
  wrapper (`formula_m = list(...)`) is **deferred** to a follow-up spec.
- **Gaussian outcome and mediators only** (medfit's parallel class has no
  `family_*` slots; non-Gaussian is blocked until medfit adds them).

---

## 2. Estimand & mathematics

For k parallel mediators, treatment X (treated `x` = `x_value`, control `x*` =
`x_ref`), contrast `δ = x_value − x_ref`. Treatment is held at `x_value` for
**both** potential outcomes; only the mediators vary between their treated
`Mⱼ(x)` and control `Mⱼ(x*)` levels (drawn independently). The direct effect
`c′` cancels.

Each mediator model: `Mⱼ = i_mⱼ + aⱼ·X + (covariates) + ε_mⱼ`,
`ε_mⱼ ~ N(0, σ²_mⱼ)`. Outcome: `Y = i_y + c′·X + Σⱼ bⱼ·Mⱼ + (covariates) + ε_y`,
`ε_y ~ N(0, σ²_Y)`.

Potential outcomes (X fixed at `x_value`):

- `Y_t = i_y + c′·x_value + Σⱼ bⱼ·Mⱼ(x)   + ε_t`
- `Y_c = i_y + c′·x_value + Σⱼ bⱼ·Mⱼ(x*)  + ε_c`

The mediator difference `Mⱼ(x) − Mⱼ(x*) ~ N(aⱼδ, 2σ²_mⱼ)`, independent across j.
Hence

```
Y_t − Y_c  ~  N( δ·Σⱼ aⱼbⱼ ,  2·Σⱼ bⱼ²σ²_mⱼ + 2σ²_Y )
```

and the **joint parallel P_med**:

```
P_med = P(Y_t > Y_c) + ½·P(Y_t = Y_c)
      = Φ(  δ · Σⱼ aⱼbⱼ  /  √( 2·Σⱼ bⱼ²σ²_mⱼ + 2σ²_Y )  )
```

This recovers the single-mediator formula `Φ(ab/√(2(b²σ²_M+σ²_Y)))` exactly at
k = 1. Reversing the contrast (`δ → −δ`) maps `P_med → 1 − P_med`.

**Indirect effect (total):** `IE = Σⱼ aⱼbⱼ` (matches medfit's
`ParallelMediationData` printed IE).

---

## 3. Input contract (medfit `ParallelMediationData`)

Verbatim S7 properties consumed (medfit `R/classes.R`):

- `a_paths` (numeric, length k): X → each Mⱼ
- `b_paths` (numeric, length k): each Mⱼ → Y (same outcome equation)
- `c_prime` (numeric scalar): X → Y direct (cancels in P_med)
- `sigma_mediators` (numeric, length k, or NULL): per-mediator residual SD
- `sigma_y` (numeric scalar, or NULL): outcome residual SD
- `estimates` (named numeric) + `vcov` (matrix): for the parametric bootstrap
- `mediators`, `mediator_predictors` (list len k), `outcome_predictors`,
  `treatment`, `outcome`, `data`, `n_obs`

**Gaussian indicator:** the parallel class carries **no** `family_m`/`family_y`.
The Gaussian guard is therefore "`sigma_mediators` and `sigma_y` are all present
and finite". If any is NA/NULL → error pointing at the medfit family-slot
enhancement as the unblock path.

---

## 4. API / entry point

New S7 method, mirroring the existing `MediationData` method:

```r
pmed(object,                  # medfit::ParallelMediationData
     x_ref = 0, x_value = 1,
     method = c("parametric_bootstrap", "nonparametric_bootstrap",
                "plugin", "mbco"),
     n_boot = 1000, ci_level = 0.95, seed = NULL, ...) -> PmedResult
```

Returns the **existing `PmedResult`** (no class change). Mapping:
`estimate` = joint P_med; `ie_estimate` = Σaⱼbⱼ; `ie_ci_*`, `boot_estimates`,
`ie_boot_estimates` as for single-mediator; `source_extract` holds the parallel
object; `method` string unchanged.

Reached via:

```r
ex <- medfit::extract_mediation(model_m, model_y = model_y, treatment = "X",
                                mediator = c("M1", "M2"), structure = "parallel")
pmed(ex, method = "mbco")
```

**Deferred:** a formula wrapper `pmed(formula_y, formula_m = list(...))` — out of
scope for this spec.

---

## 5. Computation (Approach A: closed-form kernel)

Single kernel, shared by all four methods:

```r
.pmed_parallel_closed(a_vec, b_vec, sigma_m_vec, sigma_y, delta) =
  pnorm( delta * sum(a_vec * b_vec) /
         sqrt( 2 * sum(b_vec^2 * sigma_m_vec^2) + 2 * sigma_y^2 ) )
```

- **plugin** — closed-form point estimate; `ie = Σaⱼbⱼ`; CI = NA.
- **parametric_bootstrap** — draw θ\* ~ N(`estimates`, `vcov`) (MASS::mvrnorm);
  per draw pull `a_vec`, `b_vec`, σ's **by name** from the design, evaluate the
  kernel and Σaⱼbⱼ; quantile CIs. Reuses the `R/compute-bootstrap.R` pattern.
- **nonparametric_bootstrap** — resample rows with replacement, re-extract the
  parallel structure (refit k mediator models + outcome model), evaluate the
  kernel per resample.
- **mbco** — deterministic LR-inversion interval (heaviest increment). Reuse the
  **generic** `.mbco_invert` engine. New parallel-specific pieces:
  - `.mbco_prep_parallel` — design matrices + free fits (a_hats, b_hats,
    Vm_hats, Vy_hat, ll_free, delta-method SD of Σaⱼbⱼ).
  - `.mbco_ll_constrained_pmed_parallel(prep, qstar, …)` — maximize the joint
    Gaussian log-likelihood subject to `P_med = Φ(qstar)`. The constraint
    couples `δ·Σaⱼbⱼ` and `2Σbⱼ²σ²_mⱼ + 2σ²_Y`, pinning σ²_Y; optimize the
    coupled core params (aⱼ, log Vmⱼ, log|bⱼ|) with covariates OLS-profiled;
    Nelder-Mead (penalty wall for infeasible σ²_Y).
  - `.mbco_ll_constrained_ie_parallel(prep, ie0)` — constraint `Σaⱼbⱼ = ie0`.
  - Multistart where the profile is multimodal (as single-mediator IE needed).

---

## 6. Error handling / guards

- **Non-Gaussian** (any `sigma_mediators[j]` or `sigma_y` NA/NULL) → informative
  error naming the medfit family-slot enhancement; suggests no bootstrap fallback
  exists yet for parallel non-Gaussian.
- `x_ref == x_value` (degenerate contrast) → error, as in single-mediator mbco.
- `mbco` additionally requires `@data` present (for constrained refits) → error
  if NULL.
- `< 2` mediators → error (use the single-mediator `MediationData` path); should
  not occur given medfit's validator, but guarded.

---

## 7. Testing strategy (expansive)

Estimand correctness:

- **Closed-form vs independent brute-force joint simulator** (draws all Mⱼ(x),
  Mⱼ(x*) and sums bⱼMⱼ; shares no code with the kernel) — tolerance ~1e-2.
- **k = 1 recovers** the single-mediator P_med exactly (to ~1e-10).
- **c′-invariance**: P_med unchanged as c′ varies (same a/b/σ/seed).
- **IE = Σaⱼbⱼ** (to ~1e-8).
- **Reversed contrast → 1 − P_med** (to ~1e-6).
- **Additivity sanity**: two equal mediators give the same P_med as one mediator
  with the appropriately combined a·b and variance terms.

Per method:

- plugin point estimate matches the closed form.
- parametric & nonparametric bootstrap: CIs bracket the estimate, ordered,
  reproducible under `seed`; IE-CI excludes 0 for a strong joint effect, includes
  0 for a near-null one.
- **mbco vs an independent full-dimensional oracle** (Approach-A-style oracle
  that optimizes the *un-profiled* joint likelihood) — estimate agreement to
  ~1e-6; reversed-contrast reflection; determinism (seed-free identical calls);
  `.mbco_invert` boundary/NA edge cases (already covered for the generic engine).

Full cycle: `devtools::document()` → `devtools::test()` → `R CMD check`
(0 errors, 0 warnings; the worktree `.git` NOTE is the known spurious one).
Lint-clean (`.lintr`); spelling clean (extend `inst/WORDLIST`).

---

## 8. Documentation deliverable

- **New article/vignette** `vignettes/parallel-mediation.qmd`:
  - Explains the parallel mediation structure and when to use it.
  - **Derives the joint estimand** (the `Y_t − Y_c` distribution → the closed
    form) with enough narrative for a methods reader.
  - Demonstrates all four methods on a worked example, including the
    `medfit::extract_mediation(..., structure = "parallel")` → `pmed()` flow.
  - Notes the Gaussian-only limitation and the path-specific / non-Gaussian
    extensions as future work.
  - Renders to the **pkgdown site** (article + `docs/dev/`).
- Update `NEWS.md` (Features), `pmed()` `@details`/help, reference index.
- Rebuild the pkgdown site after implementation.

---

## 9. Increment ordering (for the plan)

1. Estimand kernel `.pmed_parallel_closed` + `plugin` + `pmed(ParallelMediationData)`
   dispatch + guards. Tests: closed-form/sim oracle, k=1, c′-invariance, IE.
2. Parametric bootstrap (by-name draws from vcov). Tests: bracket/order/seed.
3. Nonparametric bootstrap (resample + re-extract). Tests: bracket/order/seed.
4. MBCO parallel (prep + constrained-pmed + constrained-ie + multistart) with the
   independent oracle. Tests: oracle agreement, reflection, determinism.
5. Vignette/article + NEWS/help/WORDLIST; rebuild site. Full-cycle check.

---

## 10. Out of scope (explicit)

- Serial/sequential P_med (separate spec).
- Path-specific per-mediator P_med (separate spec).
- Non-Gaussian parallel mediation (blocked on medfit family-slot enhancement).
- Formula convenience wrapper for parallel (`formula_m = list(...)`).
- Exposure-mediator interactions in the parallel structure.
