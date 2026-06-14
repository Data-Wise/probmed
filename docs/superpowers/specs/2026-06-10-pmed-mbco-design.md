# Design: `method = "mbco"` for `pmed()`

**Date:** 2026-06-10
**Status:** Approved (2026-06-10) — ready for implementation plan
**Package:** probmed (Data-Wise/probmed)
**Roadmap item:** "`method = "mbco"` — closed-form MBCO interval (prototype in `research/pmed/code/mbco_pmed.R`)"

---

## 1. Purpose

Add a third inference method to `pmed()` that produces a confidence interval for
P_med (and the indirect effect `a·b`) by **likelihood-ratio inversion** — the
Model-Based Constrained Optimization (MBCO) procedure of Tofighi & Kelley (2020),
already prototyped for the continuous single-mediator model in
`~/projects/research/pmed/code/mbco_pmed.R`.

Unlike the existing `parametric_bootstrap` / `nonparametric_bootstrap` methods,
MBCO is **not a resampling method**: the interval is the set of P_med values that
a Wilks likelihood-ratio test does not reject, located by root-finding. It is
therefore deterministic (seed-free) and grid-resolution-independent.

## 2. Statistical background

Gaussian single-mediator model with optional covariates `C`:

```
M = iM + a·X + g'C + eM,        Var(eM) = Vm
Y = iY + cp·X + b·M + h'C + eY, Var(eY) = Vy
```

P_med reduces to a closed form in this model:

```
P_med = Φ(d),   d = a·δ·b / sqrt(2·(b²·Vm + Vy)),   δ = x_value − x_ref
```

The covariate coefficients `g, h` and the direct effect `cp` cancel in the
potential-outcome contrast (the mediator is drawn at treated and control levels
at the *same* covariate values, treatment held fixed), so they are **nuisance
parameters** — they do not enter the P_med formula but they do change the MLEs of
`a, b, Vm, Vy`.

**MBCO constraint (P_med).** Imposing `P_med = p*` is equivalent to `d = q* =
qnorm(p*)`. Solving for the outcome variance:

```
Vy = (a·δ·b)² / (2·q*²) − b²·Vm        (q* ≠ 0)
```

Fixing `Vy` by this formula turns the nonlinear constraint into an ordinary
optimization. The MBCO profile log-likelihood at `p*` is the maximized joint
likelihood subject to this pin. The interval is `{ p* : −2(L0(p*) − L1) ≤
χ²_{1,level} }`, with endpoints found by `uniroot`.

**MBCO constraint (indirect effect).** Structurally identical: impose `a·b = ie₀`
(pin `b = ie₀ / a`, sign-managed; `Vy` free), maximize, invert the same Wilks
test. This is the MBCO test for the mediated effect — the same procedure
`RMediation::mbco` performs, reimplemented inline for our specific Gaussian
linear model to keep probmed self-contained.

## 3. Scope

**Supported**

- Gaussian outcome **and** Gaussian mediator.
- Single mediator.
- Covariates in either or both equations.
- Arbitrary treatment contrast `x_ref` / `x_value` (via `δ = x_value − x_ref`).

**Out of scope (YAGNI)**

- Non-Gaussian / binary families — the closed-form `Vy` solve is Gaussian-only;
  the bootstrap methods remain the route for those.
- Multi-mediator (parallel / serial) structures — separate `MediationData`
  subclasses that `pmed()` does not yet handle.
- A lavaan dependency / delegation to `RMediation::mbco` (its API takes fitted
  lavaan nested models, which would pull lavaan into the critical path).

**Guards (clear errors)**

- `family_y` or `family_m` non-Gaussian → error naming the family and pointing to
  `method = "parametric_bootstrap"` / `"nonparametric_bootstrap"`.
- `extract@data` is `NULL` → error explaining MBCO refits constrained models and
  needs the raw data (available when `pmed()` is called via the formula
  interface, or when the `MediationData` carries `@data`).

## 4. Architecture

New file `R/compute-mbco.R`.

### 4.1 One inversion engine, two estimands

```
.mbco_interval(ll_free, ll_constrained, est, level)
```

Performs the LR inversion (step outward from `est` until the statistic exceeds
the χ² cutoff, then `uniroot` the crossing), carried over from the prototype's
`endpoint()` / `excess()` logic. Resolution-independent; returns `c(lower,
upper)`.

Two callers supply different `ll_constrained` closures:

- **P_med:** constraint `d = qnorm(p*)`, `Vy` pinned by the formula in §2.
- **IE:** constraint `a·b = ie₀`, `b` pinned to `ie₀ / a`, `Vy` free.

Both reuse `.ll_free` — the sum of the two unconstrained Gaussian regression
log-likelihoods (`M ~ X + C`, `Y ~ X + M + C`).

### 4.2 Constrained fit with covariates — profiled likelihood (Approach B)

Numerically optimize **only the small set of coupled core parameters** —
fixed-dimension regardless of covariate count. The exact core depends on which
constraint is active:

- **P_med constraint:** optimize `(a, b, Vm)` (3-D); `Vy` pinned by the formula.
- **IE constraint:** optimize `(a, Vm)` with `b = ie₀ / a` pinned; `Vy`
  free-MLE'd from the Y-equation residuals.

Given the core values:

- the M-equation intercept and covariate coefficients are closed-form OLS of the
  partial residual `M − a·X` on `[1, C]`;
- the Y-equation intercept, `cp`, and covariate coefficients are closed-form OLS
  of `Y − b·M` on `[1, X, C]` (GLS reduces to OLS with a scalar `Vy`);
- `Vy` is pinned by the active constraint (P_med) or free-MLE'd (IE).

Design matrices are built from `extract@mediator_predictors` and
`extract@outcome_predictors` against `extract@data`.

**Optimizer.** `optim` **Nelder-Mead**, warm-started from the free OLS/MLE fits.
Nelder-Mead is chosen deliberately over a gradient method: the constrained
objective has a discontinuous penalty wall where `Vy ≤ 0` becomes infeasible
(§4.2 edge cases), which breaks BFGS gradient/line-search steps — this is why the
prototype uses Nelder-Mead. Parameterization keeps the search well-behaved:
optimize `log(Vm)` (positivity) and a sign-pinned `b = sign(target)·exp(·)` so
the fit targets `p*` rather than `1 − p*`.

**Edge cases (from prototype).**

- `q* ≈ 0` (P_med = 0.5): null indirect submodel with `b = 0` (`Y ~ X + C`).
- `Vy ≤ 0`: infeasible → treated as rejected (push the LR statistic above the
  cutoff).
- `sign(a·b)` pinned to the sign of the target so the profile tracks the correct
  branch.

### 4.3 Output (`PmedResult`)

| Property | Value |
|---|---|
| `method` | `"mbco"` |
| `estimate` | `Φ(d̂)` at the free MLEs |
| `ci_lower` / `ci_upper` | P_med interval endpoints |
| `ci_level` | `level` |
| `ie_estimate` | `a·b` |
| `ie_ci_lower` / `ie_ci_upper` | IE interval endpoints |
| `n_boot` | `NA_integer_` |
| `boot_estimates` / `ie_boot_estimates` | `numeric(0)` |
| `converged` | optimizer success flag |
| `x_ref` / `x_value`, `source_extract` | as for other methods |

The existing print/summary methods already render `NA` `n_boot` and an empty
bootstrap vector (the `plugin` method exercises this path today).

### 4.4 Wiring

- `R/methods-pmed.R`: add `"mbco"` to the `method` `match.arg(...)` vector in
  **both** `pmed()` methods (formula and `MediationData`).
- `R/compute-core.R`: add an `mbco =` branch to the `switch()` in
  `.pmed_compute`, dispatching to `.pmed_mbco(extract, x_ref, x_value, ci_level,
  ...)`.

## 5. Testing — `tests/testthat/test-pmed-mbco.R`

1. **No-covariate agreement.** Point estimate and interval match the standalone
   prototype within tolerance. (No `RMediation::mbco` test: its API takes fitted
   lavaan `h0`/`h1` models, so an automated cross-check would couple the suite to
   lavaan + RMediation for marginal gain — the prototype and the §5.2 brute-force
   grid are sufficient oracles. A manual `RMediation` sanity check can be run
   out-of-suite during development.)
2. **Covariate correctness.** Profiled-likelihood interval matches a brute-force
   grid LR inversion on a model with one or more covariates.
3. **Coverage sanity.** Interval brackets the estimate; endpoints in (0, 1);
   `ci_lower ≤ estimate ≤ ci_upper`.
4. **IE interval.** Brackets `a·b`; excludes 0 when the bootstrap IE CI does;
   includes 0 at a near-null model.
5. **Guards.** Binary `family_y` errors with an informative message; `NULL`
   `@data` errors with an informative message.
6. **Determinism.** MBCO is seed-free — two calls on the same input return
   identical endpoints.

## 6. Documentation

- `R/methods-pmed.R`: extend the `@param method` docs to list `"mbco"`.
- `CLAUDE.md`: move MBCO from the "Future" roadmap list to "Completed"; add a row
  to the Computation Methods table.
- `vignettes/introduction.qmd`: add an MBCO subsection alongside the existing
  `method = "parametric_bootstrap"` discussion (this is the vignette that teaches
  the inference-method choices).

## 7. Dependencies

**No new dependency of any kind.** `stats` (already imported) covers `optim`,
`lm`/OLS, `qnorm`, `qchisq`, `uniroot`, `pnorm`. `RMediation` is **not** added —
neither Imports nor Suggests — for the reason in §5.1.

## 8. Files touched

| File | Change |
|---|---|
| `R/compute-mbco.R` | **new** — engine + `.pmed_mbco()` + constrained fits |
| `R/compute-core.R` | add `mbco` branch to `.pmed_compute` switch |
| `R/methods-pmed.R` | add `"mbco"` to `match.arg`; doc the param |
| `tests/testthat/test-pmed-mbco.R` | **new** — tests above |
| `CLAUDE.md` | roadmap + method table update |
| `vignettes/introduction.qmd` | MBCO subsection by the method discussion |
