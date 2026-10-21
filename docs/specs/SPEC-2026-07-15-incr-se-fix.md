# SPEC — Fix `incr_pmed()` g-score SE mis-weighting (Issue #20)

- **Date:** 2026-07-15
- **Issue:** [#20](https://github.com/Data-Wise/probmed/issues/20)
- **Type:** Bug fix (variance-scale) + docstring sync + regression test
- **Target branch:** `feature/incr-se-fix` → PR to `dev`
- **Estimand impact:** none (point estimate unchanged — T2 is mean-zero)
- **Scope:** one-line code change + Roxygen formula + one new test

---

## Problem

`incr_pmed()` standard errors over-cover, worsening with δ (empirical se_ratio
1.37 → 4.40 across δ = 0.5 → 2.0). Root cause: the efficient-IF g-score term
uses the **bare** corner contrasts `gamma_med`/`gamma_dir`, while the point
estimate and base IF term carry the `qp = dq/dδ` weight. Since `tot` shrinks
with δ but the bare-γ numerator does not, the term inflates ~δ/(g(1−g)).

Confirmed against source at [R/incremental-pmed.R:139](../../R/incremental-pmed.R):

```r
# current (buggy) — bare gamma
psi_gscore <- dqg * resid * (dir * gamma_med - med * gamma_dir) / tot^2
```

`a_med`/`a_dir` (the qp-weighted contrasts) are already computed at lines
131–132, so the fix is local and free.

## Fix

```r
# fixed — qp-weighted, matches point estimate + psi_base
psi_gscore <- dqg * resid * (dir * a_med - med * a_dir) / tot^2
```

Issue-reported calibration after fix: se_ratio 1.006 / 1.081 / 1.147 at
δ = 0.5 / 1.0 / 2.0 (vs 1.35 / 2.30 / 4.57 before).

## Tasks

| # | Task | File | Acceptance |
|---|------|------|------------|
| 1 | Swap bare γ → `a_med`/`a_dir` in `psi_gscore` | `R/incremental-pmed.R:139` | line matches fix above; `devtools::load_all()` clean |
| 2 | Update Roxygen g-score formula to qp-weighted form | `R/incremental-pmed.R:67` | docstring no longer states bare `gamma_med - gamma_dir`; `devtools::document()` regenerates Rd with no other diff |
| 3 | Add SE-calibration regression test | `tests/testthat/test-incremental-pmed.R` | `se_ratio ∈ [0.9, 1.20]` on the issue DGM at δ = c(0.5, 1, 2); fails on pre-fix code |
| 4 | NEWS bullet under a `# probmed 0.3.0.9000` (dev) heading | `NEWS.md` | one line crediting the SE fix + issue #20 |
| 5 | Full `R CMD check --as-cran` in the worktree | — | 0E / 0W / 1 NOTE (existing Remotes-pin NOTE only) |

### Task 3 — test design notes

- Reuse the issue's synthetic DGM (n = 1000, τ = 0). Keep reps modest
  (≈150–200) to bound runtime; wrap in `skip_on_cran()` if wall-time > a few s.
- Assert `se_ratio = mean(se_hat) / empSD(Pmed_hat)` within `[0.9, 1.20]`.
  The upper bound is 1.20 (not 1.15) to absorb the known residual ~1.147 at
  δ = 2 plus Monte-Carlo noise at low reps — tighten only if reps are raised.
- **Anti-regression check:** confirm the test FAILS against the pre-fix line
  (se_ratio > 2 at δ = 2) before committing the fix — a test that can't fail
  on the bug is decoration (e2e-before-pr).

## Verification / evidence to record in PR body

- `devtools::test()` — full pass count, new test green.
- The pre-fix se_ratio table (bug reproduced) vs post-fix table (calibrated),
  quoted from the actual test run, not asserted.
- `R CMD check --as-cran` counts.

## Out of scope

- The residual ~1.15 at δ = 2 (defensible band; a further scaling refinement is
  a separate investigation, not this fix).
- Any change to `ward_residual()` / gauge or Wasserstein surfaces.

## Rollback

Single-line revert of Task 1 + docstring; no data or API surface changes.
