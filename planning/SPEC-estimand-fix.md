# SPEC: Fix probmed's core P_med estimand to match the mediation definition

> **Status:** IMPLEMENTED on `feature/estimand-fix` (medfit family slot on
> `feature/family-slot`, v0.2.0.9000). Gaussian core matches the closed form to
> MC error, c' invariance holds, binary outcomes are non-degenerate, and both
> bootstraps were repaired (parametric was mis-indexing intercepts as a/b).
> The pmed manuscript's §6 illustrative code depends on this fix (pointer kept in
> `pmed/docs/SPEC-probmed-estimand-fix.md`).

## 1. Summary

`probmed::pmed()` does **not** compute the estimand the manuscript defines.
Its core routine `.pmed_core_simple()` holds the **mediator** fixed and varies
**treatment X** — a direct-effect-style contrast whose value depends on the
direct effect `c'`. The manuscript's P_med holds **X fixed at treatment** and
varies the **mediator** between its treated and control values — the mediation
(indirect) contrast, independent of `c'`. The two give different numbers, often on
opposite sides of 0.5.

## 2. The two estimands

**Manuscript (Definition 1; implemented in `research/pmed/code/simulation_functions.R::pmed_point` and `code/mbco_pmed.R`):**
$$P_{med} = P\big(Y(1, M(1)) > Y(1, M(0))\big) + \tfrac12 P\big(Y(1,M(1)) = Y(1,M(0))\big).$$
Fix $X=1$; draw $M(1)\sim$ treated and $M(0)\sim$ control mediator distributions
**independently**; both outcomes evaluated at $X=1$. The direct effect $c'$ enters
both outcomes identically and **cancels** — P_med is a function of $a, b, \sigma_M, \sigma_Y$ only.

**probmed (current code + `R/generics.R` + `R/aaa-imports.R` + `CLAUDE.md`):**
$$P\big(Y_{x^*, M_x} > Y_{x, M_x}\big) = P\big(Y(0, M(1)) > Y(1, M(1))\big).$$
Draw $M$ under **treatment only**; vary $X$ between $x^*$ and $x$. Here $c'$ does
**not** cancel, so the value is essentially $\Phi(-c'/\sqrt2)$ — a direct-effect quantity.

## 3. Evidence

- `R/compute-core.R` `.pmed_core_simple()`: `m_x <- rnorm(n_sim, a*x_value, sigma_m)`
  is drawn **once** (treated); `mu_m_xref <- a*x_ref` is **computed but never used**;
  the two outcomes differ only in `c'*x_ref` vs `c'*x_value` (X varies, M held).
- Both bootstraps (`.pmed_parametric_boot`, `.pmed_nonparametric_boot`) call
  `.pmed_core_simple()`, so they inherit the same contrast.
- **Concrete diff** (a=0.5, b=0.5, c'=0.3, n=2e5): manuscript P_med = **0.564**
  (analytic $\Phi(ab/\sqrt{2(b^2+1)})$ = 0.563); probmed contrast = **0.416**
  (analytic $\Phi(-c'/\sqrt2)$). Different magnitude, opposite side of 0.5.
- Diagnostic: probmed's value changes with `c'` and is unchanged by swapping which
  mediator level you call "treated"; a true mediation P_med is the opposite.

## 4. Root cause & the correct algorithm

`.pmed_core_simple()` must draw **both** mediator levels independently and hold X at
the treatment value for both outcomes:

```r
# corrected core (Gaussian outcome):
mu_m_t <- a * x_value           # E[M | X = x_value]  (mediator intercept cancels in the contrast)
mu_m_c <- a * x_ref             # E[M | X = x_ref]
m_t <- rnorm(n_sim, mu_m_t, sigma_m)        # M(1)  (independent draws)
m_c <- rnorm(n_sim, mu_m_c, sigma_m)        # M(0)
# BOTH outcomes at X = x_value  -> the direct-effect term c'*x_value cancels
y_t <- b * m_t + c_prime * x_value + rnorm(n_sim, 0, sigma_y)   # Y(1, M(1))
y_c <- b * m_c + c_prime * x_value + rnorm(n_sim, 0, sigma_y)   # Y(1, M(0))
pmed <- mean(y_t > y_c) + 0.5 * mean(y_t == y_c)               # tie term (Definition 1)
```
Binary/non-Gaussian Y (the `else` branch): same change — draw both mediator levels,
form the linear predictors at `X = x_value`, map through the link / draw Bernoulli
responses, and **include the 0.5 tie term** (currently omitted).

## 5. Affected files (in probmed)

| File | Change |
|------|--------|
| `R/compute-core.R` | rewrite `.pmed_core_simple()` (vary M, fix X=x_value, add tie term) |
| `R/compute-bootstrap.R` | no logic change (calls the fixed core) — re-verify |
| `R/generics.R`, `R/aaa-imports.R`, `R/probmed-package.R`, `CLAUDE.md` | correct the documented estimand to $P(Y(x,M(x)) > Y(x,M(x^*)))$ |
| `tests/testthat/` | add a test pinning P_med to the mediation value (see §6); existing IE tests unaffected (a*b unchanged) |

## 6. Validation

- **Agreement:** `probmed::pmed(method="plugin")` must match `research/pmed/code/.../pmed_point()`
  and the closed form $\Phi(ab/\sqrt{2(b^2\sigma_M^2+\sigma_Y^2)})$ to MC error on a
  Gaussian DGM (target 0.564 for a=b=.5, c'=.3).
- **No c' dependence:** P_med must be invariant to changing `c'` (holding a, b, σ).
- **Manuscript example:** on `memory_exp`, `probmed::pmed()` must reproduce
  $\hat P_{med}=0.68$ (currently only the self-contained code does).
- `devtools::test()` + `R CMD check`.

## 7. Downstream / sequencing

1. **First** fix the core estimand (this spec) on a probmed feature branch → PR to `dev`.
2. **Then** the manuscript §6 code chunk (`probmed::pmed(...)`) will actually
   reproduce the reported 0.68; no manuscript prose change needed (the reported
   numbers already use the correct estimand via the self-contained code).
3. **Then** port the closed-form MBCO into probmed as `method="mbco"`. Prototype:
   `research/pmed/code/mbco_pmed.R`.

## 8. Notes
- probmed did not load on the manuscript author's Mac (Homebrew R compiles deps from
  source; `medfit` absent) — analysis used probmed's exact analytic contrast. Re-run
  the §6 diff on a machine with `medfit` installed before/after the fix.
- This is a behavior-changing fix to a pre-release package (v0.0.0.9000); no CRAN
  users affected.
