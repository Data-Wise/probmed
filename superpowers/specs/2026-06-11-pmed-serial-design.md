# Spec: Serial (joint) P_med — `method` support for `SerialMediationData`

**Date:** 2026-06-11
**Status:** approved (brainstorm) → ready for plan
**Package:** probmed
**Depends on:** medfit (`SerialMediationData`, Gaussian extraction)
**Precedent:** parallel joint P_med (PR #3) — reuse its kernel/method structure.

---

## 1. Summary

Extend P_med to a **serial (sequential) mediator chain**
`X → M₁ → M₂ → … → M_k → Y` (a *pure* chain: each Mⱼ depends only on its
predecessor; Y depends only on M_k). This spec covers the **joint** estimand:
P_med for moving the whole chain from its control to its treated realization.

Entry via a new S7 method `pmed(<medfit::SerialMediationData>, …)` returning the
existing `PmedResult`. All four methods: `plugin`, `parametric_bootstrap`,
`nonparametric_bootstrap`, `mbco`.

**Scope decisions (locked in brainstorm):**

- **Serial only**, **joint estimand only** (path-specific deferred).
- **All four methods.**
- **Closed-form kernel** (Gaussian exact; no simulation), as in parallel.
- **Object-entry only**; formula wrapper deferred.
- **Gaussian outcome and mediators only** (medfit's serial class has no
  `family_*` slots; non-Gaussian blocked until medfit adds them).
- **Shared MBCO helpers:** rename `.mbco_parallel_starts` / `_optim_best` →
  `.mbco_starts` / `.mbco_optim_best` (structure-agnostic) so serial and
  parallel share one multistart/optim core. `.mbco_invert`, `.qr_sse`,
  `.mbco_gll`, `.ols_sse` are already generic and reused as-is.

---

## 2. Estimand & mathematics

Pure chain with Gaussian equations:

- `M₁ = i₁ + a·X + ε₁`, `ε₁ ~ N(0, V₁)`
- `Mⱼ = iⱼ + d_{j−1}·M_{j−1} + εⱼ`, `εⱼ ~ N(0, Vⱼ)` (j = 2..k)
- `Y = i_y + c′·X + b·M_k + ε_Y`, `ε_Y ~ N(0, V_Y)`

Treatment is held at `x = x_value` for **both** potential outcomes; only the
mediator chain varies between its treated realization (generated with X = x at
the M₁ stage) and its control realization (X = x* at the M₁ stage), drawn
independently. The direct effect `c′` cancels. With `δ = x_value − x_ref`, the
chain difference propagates:

```
ΔM₁ ~ N(a·δ, 2V₁)
ΔMⱼ = d_{j−1}·ΔM_{j−1} + Δεⱼ      Var[ΔMⱼ] = d_{j−1}²·Var[ΔM_{j−1}] + 2Vⱼ
ΔY  = b·ΔM_k + Δε_Y
```

so `ΔY ~ N( a·(∏d)·b·δ,  2·(b²·Var[ΔM_k] + V_Y) )` where `Var[ΔM_k]` is the
recursion above (closed under the chain). The **joint serial P_med**:

```
P_med = Φ(  a·(∏d)·b·δ  /  √( 2·(b²·Var[ΔM_k] + V_Y) )  )
```

Recovers the single-mediator formula at k = 1 (no `d`, `Var[ΔM₁] = 2V₁` →
`Φ(ab δ / √(2(b²V₁ + V_Y)))`). Reversing the contrast maps `P_med → 1 − P_med`.

**Indirect effect:** `IE = a·(∏d)·b` (matches medfit's serial IE).

Helper: a small recursive function `var_delta_chain(d_vec, V_m_vec)` returns
`Var[ΔM_k]`; the kernel is then a one-liner.

---

## 3. Input contract (medfit `SerialMediationData`)

Verbatim S7 properties consumed:

- `a_path` (numeric scalar): X → M₁
- `d_path` (numeric, length k−1): the chain links M₁→M₂, …, M_{k−1}→M_k
- `b_path` (numeric scalar): M_k → Y
- `c_prime` (numeric scalar): X → Y (cancels)
- `sigma_mediators` (numeric, length k, or NULL): per-mediator residual SD
- `sigma_y` (numeric scalar, or NULL)
- `estimates` (named numeric) + `vcov` (matrix) for the parametric bootstrap
- `mediators` (length k, chain order), `mediator_predictors` (list len k),
  `outcome_predictors`, `treatment`, `outcome`, `data`, `n_obs`

**Gaussian indicator:** no `family_*` slots → require all `sigma_mediators` and
`sigma_y` present and finite, else error (pointing at the medfit family-slot
enhancement). Same heuristic and lavaan caveat as parallel.

---

## 4. API / entry point

```r
pmed(object,                  # medfit::SerialMediationData
     x_ref = 0, x_value = 1,
     method = c("parametric_bootstrap", "nonparametric_bootstrap",
                "plugin", "mbco"),
     n_boot = 1000, ci_level = 0.95, seed = NULL, ...) -> PmedResult
```

Returns the existing `PmedResult` (`estimate` = joint P_med; `ie_estimate` =
`a·∏d·b`; CI/IE-CI/boot vectors as in parallel; `source_extract` holds the
serial object). Registered in `R/methods-pmed.R` (after the `pmed` generic, as
the parallel method is). Reached via
`medfit::extract_mediation(..., structure = "serial")`.

Formula wrapper deferred.

---

## 5. Computation (closed-form kernel + four methods)

- **Kernel** `.pmed_serial_closed(a, d_vec, b, V_m_vec, V_y, delta)` →
  `pnorm(a*prod(d)*b*delta / sqrt(2*(b^2*var_delta_chain(d_vec, V_m_vec) + V_y)))`.
- **plugin** — closed-form point estimate; `ie = a*prod(d)*b`; CI = NA.
- **parametric_bootstrap** — draw the structural coefficients from
  `N(estimates, vcov)` using the **source-agnostic structural aliases** the
  medfit serial extractor exposes: `a`, `d1`…`d_{k−1}`, `b` (NOT lm-prefixed
  `m{j}_…` names — the parallel review showed those break on lavaan extracts).
  Guard on any missing alias name. Residual SDs held fixed. Evaluate the kernel
  + IE per draw; quantile CIs.
- **nonparametric_bootstrap** — resample rows, refit the chain equations
  (`M₁ ~ X(+cov)`, `Mⱼ ~ M_{j−1}(+cov)`, `Y ~ X + M_k(+cov)`) from
  `mediator_predictors`/`outcome_predictors`, evaluate the kernel per resample.
- **mbco** — deterministic LR inversion (heaviest increment). Reuse `.mbco_invert`.
  New: `.mbco_prep_serial` (chain design matrices + free fits + QR precompute +
  delta-method SD of `a·∏d·b`); `.mbco_ll_constrained_pmed_serial` (constraint
  pins `V_Y` via the P_med formula; optimize `a, d_1..d_{k−1}, b, log V_j` with
  intercepts/covariates OLS-profiled through precomputed QR; Nelder-Mead with a
  penalty wall on infeasible `V_Y`; deterministic multistart via `.mbco_starts`);
  `.mbco_ll_constrained_ie_serial` (constraint `a·∏d·b = ie0`, pin one path
  coord, e.g. `b = ie0 / (a·∏d)`, with multistart). Performance carries the
  parallel lessons: QR-precomputed `.qr_sse`, `reltol = 1e-8`, coarse IE
  inversion step.

---

## 6. Error handling / guards

- Non-Gaussian (any `sigma_mediators[j]`/`sigma_y` NA/empty) → informative error.
- `x_ref == x_value` → error (mbco).
- `mbco`/nonparametric require `@data`.
- `< 2` mediators → error (use single-mediator `MediationData`).
- Degenerate chain (`a`, any `d_j`, or `b` ≈ 0) → IE/MBCO interval may be NA;
  reported separately, documented (as parallel does for `converged`).

---

## 7. Testing strategy (expansive)

- Closed-form vs an **independent brute-force chain simulator** (draws the full
  treated/control chains, no shared code with the kernel) — tolerance ~1e-2.
- **k = 1 recovers** single-mediator P_med exactly; **c′-invariance**;
  **IE = a·∏d·b**; **reversed contrast → 1 − P_med**.
- Per method: bootstrap CIs bracket/order/seed-reproducible; IE-CI excludes 0
  for a strong chain, includes 0 near null.
- **mbco vs an independent full-dimensional oracle** (un-profiled chain
  likelihood) at **k = 2 and k = 3**; determinism; constrained-LL-at-MLE ≈ free.
- Guards: non-Gaussian (mediator-side + non-plugin entry), degenerate contrast,
  missing data; lavaan-sourced parametric bootstrap (skip if lavaan absent).
- Full cycle: `document → test → check` (0 errors/0 warnings; spurious `.git`
  NOTE only); lint + spelling clean.

---

## 8. Documentation deliverable

- **New vignette** `vignettes/serial-mediation.qmd`: derives the recursive chain
  variance and the closed form, demonstrates all four methods on a worked
  `extract_mediation(..., structure = "serial")` example, notes the Gaussian /
  joint-only / pure-chain limits.
- **Register the article in `_pkgdown.yml`** — both the navbar menu and the
  `articles:` index (a missing entry fails `build_articles_index`, as the
  parallel build surfaced).
- Update `NEWS.md` (Features), `pmed()` help, `inst/WORDLIST`.
- **Install the updated package before building the site** (the standalone
  vignette render uses the installed probmed), then rebuild the pkgdown site.

---

## 9. Increment ordering (for the plan)

0. Rename `.mbco_parallel_starts`/`_optim_best` → `.mbco_starts`/`_optim_best`
   (update parallel call sites; tests stay green) — shared core for serial.
1. Kernel `.pmed_serial_closed` + `var_delta_chain` + `plugin` +
   `pmed(SerialMediationData)` dispatch + guards. Tests: simulator, k=1, c′, IE.
2. Parametric bootstrap (alias-based draws + guard). Tests: bracket/order/seed.
3. Nonparametric bootstrap (resample + refit chain). Tests: bracket/order/seed.
4. MBCO serial (prep + constrained-pmed + constrained-ie + multistart) with the
   independent oracle (k=2, k=3). Tests: oracle agreement, reflection, determinism.
5. Vignette/article + `_pkgdown.yml` registration + NEWS/help/WORDLIST; install +
   rebuild site. Full-cycle check.

---

## 10. Out of scope (explicit)

- Path-specific serial P_med (separate spec).
- Non-Gaussian serial mediation (blocked on medfit family-slot enhancement).
- Skip/direct paths beyond the pure chain (X→Mⱼ for j>1, Mᵢ→Y for i<k).
- Formula convenience wrapper.
- Parallel ⇄ serial unification beyond the shared MBCO helpers.
