# probmed repair: task list

Plan: [`tasks/plan.md`](plan.md). Technical requirements: the Codex
`planning/repair/PLAN-…` and `SPEC-…` (2026-10-09) on branch `codex/probmed-repair-2026-10-09`.
Focused tests: `devtools::test(filter = "<name>")`. Full suite: `devtools::test()`.
Check: `rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"), error_on = "warning")`.
Quick local gate: `scripts/check-local.sh`.

---

## Status (2026-10-14)

| Phase | Tasks | Where | State |
|---|---|---|---|
| 0–1 | T1–T5 | PR #41 → dev `1c39d15` | merged |
| 2 | T6–T8 | PR #42 → dev `22edf49` | merged |
| 3 | T9–T11 | PR #43 → dev `7847db5` | merged |
| 4 | T12–T14 | PR #44 → dev `0fb9e81` | merged |
| 5 | T15–T17 | PR #45 `feature/repair-calibration` → dev | open; 17 PASS / 15 FAIL / 4 INCONCLUSIVE (`docs/reports/CALIBRATION-2026-10-14.md`) |
| side | T18 | dev `ecffb97` | done; #38 comment posted |

Checkpoints B and C passed their test/check gates but **not** the user review they call for.
Stacked PRs: after each merge, rebase the next with
`git rebase --onto dev <merged-branch> <next-branch>` and retarget it to `dev`.

## Phase 0: Secure the work

### T1: Back up the Codex branch to origin  *(needs the user's yes)*
**Description:** The two repair commits exist on one disk only. Push them so they can't be lost.
**Acceptance criteria:**
- [x] `git branch -r --contains a1e3737` lists an `origin/…` ref
**Verification:**
- [x] `git ls-remote origin codex/probmed-repair-2026-10-09` returns `a1e3737`
**Dependencies:** None · **Files:** none · **Scope:** XS

---

## Phase 1: Land batches 1–2 (target branch `feature/mbco-repair`)

### T2: Clean up and rebase
**Description:** Rebase the two commits onto current dev under a `feature/*` name. Drop the files that belong to the manuscript project and the branch's `.STATUS` edit.
**Acceptance criteria:**
- [x] Branch `feature/mbco-repair` = dev + batches 1–2; no conflicts
- [x] `planning/repair/*` and `probmed-repair-handoff-*.json` removed from the package branch (or moved per open question Q2)
- [x] `.STATUS` unchanged relative to dev
**Verification:**
- [x] `git diff dev...feature/mbco-repair --stat` shows only R/, man/, tests/, NEWS, README, DESCRIPTION, vignettes, docs/reports
- [x] Full suite: 621+ passed, 0 failed
**Dependencies:** T1 · **Files:** git-only, plus removal of 3 planning files · **Scope:** S

### T3: Measure the MBCO NA rate and runtime  *(decision point)*
**Description:** The new failure rule turns one failed grid point (in either sign case) into a missing interval. Measure how often that happens and what it costs in time, before the PR.
**Acceptance criteria:**
- [x] Script `inst/sim/mbco_na_rate.R`: 200 seeded Gaussian fixtures × {null a=0, null b=0, weak, strong} × n ∈ {100, 500}, with and without one covariate
- [x] Reports per cell: share of `status` values ok/truncated/disconnected/failed for P_med and IE, plus median and 95th-percentile wall time
- [x] Result committed as `inst/sim/results/mbco_na_rate.csv` plus a short FINDINGS note
- [x] **Decision:** if NA rate > threshold (open question Q4), do T3b before T4
**Verification:**
- [x] Positive control: a planted always-failing profile yields `failed` in 100% of runs
**Dependencies:** T2 · **Files:** 1 script, 1 CSV, 1 note · **Scope:** S

### T3b: *(conditional)* Make wrong-sign-basin fitting more robust
**Description:** If T3 shows too many failures, add starts or a tolerant rule. For example, accept a P_med profile when the failing sign case is provably dominated (its likelihood bound is below the other case's maximum).
**Acceptance criteria:**
- [x] NA rate within threshold on the T3 grid
- [x] No fixture that previously failed is now reported as `ok` with a wrong endpoint (independent nloptr oracle within 1e-4)
**Verification:**
- [x] `devtools::test(filter = "repair-profile")` green; T3 script re-run
**Dependencies:** T3 · **Files:** `R/compute-mbco-verified.R`, `tests/testthat/test-repair-profile.R` · **Scope:** S

### T4: R CMD check plus the roxygen toolchain gate
**Description:** Codex regenerated `man/` with roxygen2 8.1.0. 8.x is known to break the docs for S7 methods (memory `roxygen-version-toolchain`).
**Acceptance criteria:**
- [x] `R CMD check --as-cran`: 0 errors, 0 warnings; NOTEs match the dev baseline
- [x] No codoc or "objects in \usage without \alias" warnings; DESCRIPTION has no `Config/roxygen2/version` line
- [x] Vignette `introduction.qmd` renders
**Verification:**
- [x] If warnings appear: `remotes::install_version("roxygen2", "7.3.3")`, `devtools::document()`, re-check
**Dependencies:** T2 (T3b if triggered) · **Files:** `man/*.Rd`, DESCRIPTION · **Scope:** S

### T5: Open a PR to dev with evidence
**Description:** Open the PR with the test counts, the R CMD check result, the T3 table and an end-to-end transcript from a fresh temporary install. Flag the ML-variance (SSE/n) point-estimate change for the author's review.
**Acceptance criteria:**
- [x] PR body: exact commands, pass/fail/skip counts, NA-rate table, E2E transcript (P_med and IE intervals with `@diagnostics` statuses)
- [x] PR body states the point-estimate convention change and its typical size on the example data
- [x] Leak scan of the diff is clean (no tokens, no plain-text destructive git commands)
**Verification:**
- [x] ~~CI green on the PR~~ — no CI runs on PRs into dev; local evidence in the PR body instead
**Dependencies:** T3, T4 · **Files:** none · **Scope:** XS

### ✅ Checkpoint A: after T1–T5
- [x] PR merged into dev (merge only on the user's say-so)
- [ ] dev full suite green
- [ ] Author accepted or rejected the SSE/n point-estimate change

---

## Phase 2: Support boundaries and stochastic controls (target `feature/repair-batch3`)

### T6: Reject unsupported family/link/structure before any random draws
**Description:** `.pmed_core_simple` draws Bernoulli outcomes for every non-Gaussian family (`R/compute-core.R:119,131`). Guard so that only gaussian/identity and binomial/logit|probit pass. Interactions, transformed terms and rank-deficient designs are rejected too.
**Acceptance criteria:**
- [x] Poisson / Gamma / quasi families → informative error naming the family, raised before any RNG use (`.Random.seed` unchanged)
- [x] Interaction term in either formula → informative error
- [x] Existing Gaussian and binary tests unchanged
**Verification:**
- [x] New `tests/testthat/test-support-guard.R` (one test per rejected structure + RNG-untouched assertion)
- [x] Full suite green
**Dependencies:** Checkpoint A · **Files:** `R/compute-core.R`, `R/methods-pmed.R`, new test · **Scope:** S

### T7: Forward `n_sim` and `seed` on every stochastic path
**Description:** The bootstrap calls `.pmed_core_simple` without `n_sim` (`R/compute-bootstrap.R:42,135`). Audit the parametric, nonparametric and parallel paths. Validate `n_sim` (a positive integer, kept separate from `n_boot`).
**Acceptance criteria:**
- [x] `n_sim` reaches the core on all paths (single + `compute-parallel.R`); invalid values error out
- [x] Same `seed` gives identical results on every method; different `n_sim` changes the Monte Carlo SD in the expected direction
- [x] `@param seed` / `@param n_sim` docs state the reset-then-advance RNG behavior
**Verification:**
- [x] Mock test: a spy on `.pmed_core_simple` records the received `n_sim` on each path
- [x] Full suite green
**Dependencies:** T6 · **Files:** `R/compute-bootstrap.R`, `R/compute-parallel.R`, `R/methods-pmed.R`, test · **Scope:** M

### T8: Covariance-extraction validation
**Description:** Compare `medfit::extract_mediation()@vcov` cross-equation blocks with an independent stacked estimating-equation (or paired-refit) reference, under medfit 0.4.0 and 0.6.0. Record when zero blocks are justified.
**Acceptance criteria:**
- [x] Reference written without calling medfit or probmed helpers
- [x] Disposition recorded (`docs/reports/COVARIANCE-VALIDATION-2026-10-14.md`; medfit 0.7.0; suite also green on CRAN medfit 0.3.2) per medfit version: agree / approximation justified / bug → upstream issue (filing needs a yes, since it's another repo)
**Verification:**
- [x] Test asserts agreement within tolerance on a Gaussian fixture where the zero block is provably correct
**Dependencies:** T7 · **Files:** test, report · **Scope:** S

### ✅ Checkpoint B: after T6–T8
- [x] Full suite plus `scripts/check-local.sh` green
- [ ] Review with the user before standardization

---

## Phase 3: Covariate standardization (target `feature/repair-standardize`)

### T9: Observed-sample standardization, Gaussian path
**Description:** Evaluate both cross-world arms at the same sampled covariate row, using the full model predictors and treatment at `x_value` in both outcome arms. `n_sim` is a total-pair budget.
**Acceptance criteria:**
- [x] Gaussian additive model: matches the closed-form P_med with covariates (direct-effect and covariate terms cancel)
- [x] `c'`-invariance test still passes
**Verification:**
- [x] `devtools::test(filter = "pmed-estimand")` green plus a new standardization test
**Dependencies:** Checkpoint B · **Files:** `R/compute-core.R`, `R/compute-bootstrap.R`, test · **Scope:** M

### T10: Standardization, binary path
**Description:** Same as T9 for a Bernoulli mediator and/or outcome, where covariates do **not** cancel.
**Acceptance criteria:**
- [x] Adjusted binary example agrees with an independent g-computation reference (written without probmed helpers) within Monte Carlo error at large `n_sim`
- [x] Dropping a relevant covariate changes the answer (a positive control that the covariates are actually used)
**Verification:**
- [x] New test with a known-answer g-computation fixture
**Dependencies:** T9 · **Files:** `R/compute-core.R`, test · **Scope:** S

### T11: Binary separation detector
**Description:** Detect complete or quasi-separation in the original fit and in each bootstrap/refit draw. Record original-fit failure separately from draw failures.
**Acceptance criteria:**
- [x] A planted separable fixture is flagged; a well-conditioned one is not
- [x] Rare-outcome N = 50 stress fixture reports conditional and unconditional diagnostics without changing the data-generating process
**Verification:**
- [x] New `test-separation.R`
**Dependencies:** T10 · **Files:** `R/utils.R` or a new helper, test · **Scope:** S

### ✅ Checkpoint C: after T9–T11
- [x] Full suite green; R CMD check clean
- [ ] Review with the user: is standardization behavior acceptable as the new default for adjusted models?

---

## Phase 4: Opt-in response bootstrap (target `feature/repair-response-boot`)

### T12: Simulate and refit
**Description:** Hold observed X and C fixed. Simulate M from its fitted model, then Y from its fitted model using the simulated M. Refit both equations on the same rows and recompute standardized P_med. Opt-in via `method = "response_bootstrap"` (name to confirm).
**Acceptance criteria:**
- [x] Residual-scale uncertainty propagates (checked through the delta-method SE, which includes Var(V̂m), Var(V̂y))
- [x] Default `method` unchanged
**Verification:**
- [x] Known-answer Gaussian fixture: bootstrap SD of the P_med estimate within 20% of a large-B reference
**Dependencies:** Checkpoint C, T8 · **Files:** `R/compute-bootstrap.R`, `R/compute-core.R`, `R/methods-pmed.R` · **Scope:** M

### T13: Completion and failure diagnostics
**Acceptance criteria:**
- [x] `@diagnostics` reports completed / failed / separated draw counts; the interval comes back `NA` when the completion share is below a documented floor
**Verification:**
- [x] Forced-failure fixture is reported, not hidden (mocked core returning NA; 98/100 → incomplete, 99/100 → ok)
**Dependencies:** T12 · **Files:** `R/compute-bootstrap.R`, `R/methods-print.R`, test · **Scope:** S

### T14: Documentation
**Acceptance criteria:**
- [x] roxygen + `man/` (roxygen 7.3.3), README.qmd → README.md, NEWS, `vignettes/introduction.qmd` describe standardization, the guard and the opt-in bootstrap; no claim of calibration
**Verification:**
- [x] R CMD check clean; `pkgdown::check_pkgdown()` OK; the vignette renders
**Dependencies:** T13 · **Files:** ≤5 docs files · **Scope:** M

### ✅ Checkpoint D: after T12–T14
- [x] Full suite plus R CMD check green; PR(s) to dev opened with evidence

---

## Phase 5: Calibration  *(compute needs the user's approval)*

### T15: Cost benchmark
**Acceptance criteria:**
- [x] Per-replication wall time for each method in the spec's calibration matrix; projected core-hours for the full grid
**Dependencies:** Checkpoint D · **Scope:** S

### T16: Hopper pilot, then the full array
**Acceptance criteria:**
- [x] ~~Hopper pilot + array~~ — run locally instead (14 cores, exported snapshot of `af243ab`, `PROBMED_LOAD=source`); hopper wrappers kept in `inst/sim/calibration/hopper/` (hopper has medfit 0.3.1 and no lpSolve, so it would need a fresh library first)
- [x] The run uses an installed probmed at a recorded commit, not the global library
**Dependencies:** T15 plus the user's approval · **Scope:** M

### T17: Collate and decide
**Acceptance criteria:**
- [x] PASS / FAIL / INCONCLUSIVE per cell under the spec rule, with undercoverage/overcoverage direction and failure strata
- [x] NEWS and docs restrict claims for any cell that is not PASS
**Dependencies:** T16 · **Scope:** S

---

## Independent track

### T18: #38: collate the post-fix gauge rerun
**Acceptance criteria:**
- [x] `ssh hopper ls ~/gauge_boot_postfix/` confirms the outputs exist (or the run is re-queued)
- [x] Collated CSV committed to `inst/sim/results/`; diff against the old Table 1 numbers
- [x] Draft comment for #38 written (`docs/reports/ISSUE-38-comment-DRAFT-2026-10-14.md`; posting needs a yes)
**Dependencies:** None · **Scope:** S

### T19: Housekeeping
- [x] On the next real `.STATUS` edit: fix the stale `next:` item (1) (0.3.0 already released) and record the repair status once
- [x] `AGENTS.md`: the user decides whether to keep `4901862` or revert it
**Dependencies:** None · **Scope:** XS
