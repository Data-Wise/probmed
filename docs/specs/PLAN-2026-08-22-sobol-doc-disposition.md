# Plan: `sobol_pmed()` documentation and boundary-message disposition (post-fix)

**Date:** 2026-08-22 (night). **Source:** `FINDINGS-2026-08-22-sobol-pmedw-postfix.md`
§1. **Shape:** docs + one `message()` string + test relaxation + script promotion;
no estimator change. **Branch:** `feature/sobol-doc-disposition` off `dev` (worktree;
code files touched, so not on `dev` directly). **Effort:** ~45 min + ~15 min checks.

**Status: EXECUTED** — PR #37 merged to `dev` as `c28c31e` (2026-08-22 21:39 MDT; 8/8 CI,
suite 504/0, check 0E/0W). All 12 surfaces landed; a fresh-context review before merge
found the retired story still quoted in two unplanned places (the original `sobol_pmed()`
NEWS feature bullet; the `reps` test comment) — fixed in `20791cf`. Lesson recorded
(memory `retire-a-story-grep-its-numbers`): grep the repo for the story's numbers, not
just the planned files. Side effect: `^\.git$` added to `.Rbuildignore` (worktree `.git`
file was shipping in the tarball → check NOTE).

## What the measurement licenses

- The analytic se of `Delta_m` is exact at every effect size (0.99, CV 0.08).
- Procedure B1 (default) covers 0.952-0.996 everywhere; the Wald interval for `P`
  fails only where the oracle's does (square map); Procedure A under-covers at the
  transition (0.90-0.92); B2 adds nothing.
- The bootstrap se is 0.84-1.05x the empirical SD (slightly *anti*-conservative at
  the transition), never better than the analytic default.
- `reps > 1` has nothing left to remove (empSD already at the oracle's).

So: stop telling users to switch to the bootstrap or to `reps > 1` near the
boundary; say Procedure B with the analytic se is the valid default.

## Surfaces (each with its replacement)

| # | file:lines | current claim | replacement |
|---|---|---|---|
| 1 | `R/sobol-pmed.R:47-78` (`?SobolPmedResult` "Coverage near the boundary (A-15)", item (2) "three pinned facts") | fold-split ~80% of `Var(Delta_m_hat)`; `Delta_m_hat` normal; residual se bias ~0.8x -> B covers ~0.85; bootstrap remedy 1.25x / 0.97 | Keep item (1) (Procedure A, Leeb-Potscher) verbatim. Replace item (2) with one paragraph: post-fix (2026-08-22, corrected corner weights) the analytic se of `Delta_m` matches its sampling SD (0.99, CV 0.08) from strong effects to the null, `Delta_m` Wald coverage 0.95-0.96, and `P`'s Procedure-B1 image covers 0.952-0.996 across six cells; the only near-null pathology is the Wald interval for `P` itself (0.74 at the boundary, the same for an oracle), which B avoids by construction. Add a **History** sentence: the three pre-fix "pinned facts" were the weight bug (`.corner_fit()`, PR #34); `reps > 1` and the bootstrap remain options, not remedies. |
| 2 | `R/sobol-pmed.R:80-86` "Practical guidance" | near the boundary pass `reps > 1` and `se_method = "bootstrap"`; "validated by a separate large simulation" | The analytic Procedure-B default is valid across the measured transition (cite the FINDINGS table). `reps > 1`: reproducibility only. Bootstrap: an option; post-fix slightly anti-conservative at the transition (0.84-0.90x), so not a remedy. Replace "validated by a separate large simulation" with the script + CSV paths. |
| 3 | `R/sobol-pmed.R:122-127` class `@param se_method`, `@param reps` | "near-boundary fold-split" | drop the near-boundary claims; one line each |
| 4 | `R/sobol-pmed.R:366-369` generic `@details` | "near the boundary use `reps > 1` and `se_method = "bootstrap"` ... analytic se anti-conservative there" | "Its validity inherits that of the input `Delta_m` CI, which the analytic se delivers (post-fix; see [SobolPmedResult])." |
| 5 | `R/sobol-pmed.R:395-400` generic `@param reps`, `@param se_method` | "removes the ~80% fold-split variance"; "bootstrap valid (conservative) near the boundary, analytic anti-conservative" | `reps`: "averaging over fold draws; reproducibility, small efficiency effect". `se_method`: "bootstrap is an option (B refits); post-fix no better than analytic, slightly narrower at the transition". |
| 6 | `R/sobol-pmed.R:327-336` boundary `message()` (procedure B branch) | "analytic se_Dm is anti-conservative here ... pass se_method = 'bootstrap' ... reps > 1" | "near the V_med = 0 boundary (split test p = ...); reporting Procedure B (image of the Delta_m CI, valid here; see ?SobolPmedResult)." Keep the bootstrap sub-branch text neutral ("se from nonparametric bootstrap"). |
| 7 | `R/sobol-pmed.R:218-223, 277-284` code comments | ~80% fold-split; ~0.8x anti-conservative; bootstrap 1.25x | one-line history notes pointing at the FINDINGS; keep the RNG-contract comment intact |
| 8 | `man/SobolPmedResult.Rd:155-190`, `man/sobol_pmed.Rd:50-60, 80-90` | mirrors of 1-5 | hand-mirror (roxygen2 8.0.0 here breaks S7 method docs; `tools::checkRd` both) |
| 9 | `tests/testthat/test-sobol-pmed.R:162-175` | "bootstrap se conservative -> `fb@se_Dm > 0.9 * fa@se_Dm`" | rename to "bootstrap se is finite, positive, and of the analytic se's order"; bound `> 0.6 *` (post-fix ratio 0.84-1.05 of empSD; 0.9 would be flaky) and `< 1.6 *`; keep the `Dm` identity assertion |
| 10 | `NEWS.md` (0.3.0.9000) | — | new bullet under the weight-fix section: "`sobol_pmed()` docs corrected: the near-null se story was the weight bug; Procedure B with the analytic se is valid across the transition (six-cell check against exact truth); boundary message no longer recommends the bootstrap" |
| 11 | `inst/sim/` | scripts only in the session scratchpad | promote `sobol_shipped_check.R`, `sobol_collate.R`, `pmedw_shipped_check.R` (path-relative: `pkgload::load_all(".")`, `inst/sim/phi_decomposition.R`, `GAUGE_SIM_OUT`); `inst/sim/results/sobol_postfix_summary.csv`; FINDINGS §1/§2 "Script:" lines updated to the promoted paths |
| 12 | `.STATUS` | — | worktree bullet; FINDINGS "Status" paragraph: done |

Not touched: the estimator, Procedure A/B code, `ci_B2`, `sobol_from_theta`,
vignettes (none mention sobol), `pmed-modern/03-sobol-pmed` (cross-repo; its
manuscript inherits the same correction — a separate handoff if wanted).

## Order

1. Worktree; promote scripts + CSV first (the scratchpad is session-scoped).
2. Roxygen edits 1-5, then 6-7, then Rd mirrors (8); `tools::checkRd`.
3. Test 9; run `test-sobol-pmed.R` (NOT_CRAN), then the full suite and
   `rcmdcheck --as-cran --no-manual`.
4. NEWS, FINDINGS paths, `.STATUS`; commit; PR with the six-cell table as the E2E
   (the behavior change is a message string — the E2E is `expect_no_warning` on
   the B path plus `expect_message(..., "Procedure B")` at a boundary seed, which
   fails on `dev`'s wording).

## Risks

- The `.sobol_fit` message is matched by no test today; adding the
  `expect_message` test pins it (good) but needs a seed that lands at the
  boundary (`test-sobol-pmed.R:104-115` already finds one via `cell_null`).
- Hand-mirrored Rd drifts easily; checkRd catches syntax, not content — diff
  the roxygen block against the Rd block before committing.
- Do not "fix" the bootstrap se (0.84x at the transition) in this PR; it is an
  option, documented as such. If it is ever made the default it needs its own
  measurement.
