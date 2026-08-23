# Plan: post-weight-fix follow-ups

**Date:** 2026-08-22 (evening). **Context:** PR #34 fixed the corner-EIF weight bug
(one constant per fold since PR #8); PR #35 re-measured the gauge bootstrap arm
and gates. Three estimators and two documents still rest on pre-fix numbers.
Ordered by what can be done from this repo without a decision from the user.

**Status (2026-08-22, night):** 1 and 2 DONE — `FINDINGS-2026-08-22-sobol-pmedw-postfix.md`
(sobol: oracle efficiency, Δ_m se exact, Procedure B valid everywhere; pmedW_dr: empSD −25%,
no bias, dr bootstrap calibrated — nothing to change); the sobol doc consequence shipped
as PR #37 (`c28c31e`, via `PLAN-2026-08-22-sobol-doc-disposition.md`). 3: handoff written
(`HANDOFF-2026-08-22-gauge-manuscript-postfix.md`), tracked as issue #38, and the Table 1
rerun is RUNNING on hopper (jobs 4311422 + 4311424, `~/gauge_boot_postfix/`; see
`inst/sim/hopper/README.md` "Post-fix rerun") — the manuscript edit itself is the author's.
4: still waiting on the hold.

## 1. `sobol_pmed()` post-fix coverage — in-repo, ~1 h

`sobol_pmed()` builds on `.corner_fit()` and ships an analytic (delta-method)
se, a bootstrap se, and two CI procedures (A: gated pre-test; B: sample-split).
None has been measured since the fix; its tests are structural.

- **Truth:** exact — `sobol_from_theta(truth(s, tau, FALSE)$theta)` from
  `inst/sim/phi_decomposition.R` (continuous Y only; `sobol_pmed` is
  `binY = FALSE`).
- **Script:** `sobol_shipped_check.R <s> <tau> <n> [nrep] [B]` (scratch, then
  `inst/sim/` via worktree): per dataset, `sobol_pmed()` with
  `se_method = "analytic"` (procedures A and B from one call) and, on a subset,
  `se_method = "bootstrap"`; record `p_med`, `se`, `ci`, `ci_A`, `ci_B1`,
  `ci_wald`, `boundary`, `vmed_split_reject`.
- **Cells:** manuscript 1-4 (n = 800 / 3000 x tau = 0 / 0.8, s = 1), the
  intermediate (s = 0.5, tau = 0.4) and near-null (s = 0.2) DGPs; 250 datasets
  for the analytic arm, 100 x B = 100 for the bootstrap arm.
- **Pass:** Wald / Procedure-B coverage within +/-0.03 of 0.95 in the regular
  cells, se/empSD median 0.9-1.0 with CV near the oracle's; Procedure A's gate
  rate and boundary behavior documented, not judged.
- **Deliverable:** `FINDINGS-2026-08-22-sobol-postfix.md` + results CSV; doc
  edits to `sobol_pmed` / `SobolPmedResult` only where the numbers change what
  the docs claim.

## 2. `pmedW_dr()` post-fix check — in-repo, ~1.5 h

`pmedW_dr()` has its own fold loop (same `pa()` bug, fixed in #34) and returns a
point estimate; `wasserstein_pmed(method = "dr")` wraps it with a refit
bootstrap (`n_boot`) for the se.

- **Truth:** simulate the corner laws `Y(a, M(a'))` directly from the DGP
  equations (2e6 draws) and compute `P_med^W` with `pmedW_md()` / the 1-d W2;
  no closed form needed.
- **Measure:** bias and empSD of the `pmedW_dr` point at n = 800 / 3000, pre-fix
  vs post-fix (pre-fix by sourcing `git show a73efc8^:R/wasserstein-pmed.R`),
  then bootstrap coverage on 100 datasets x `n_boot = 100` (cost is the
  constraint: ~1 s per refit).
- **Pass:** bias within MC noise; post-fix empSD <= pre-fix; bootstrap coverage
  within +/-0.04 of 0.95 away from the `NIE^W = 0` boundary.
- **Deliverable:** section in the same FINDINGS doc; `wasserstein_pmed` docs
  touched only if the bootstrap arm misbehaves.

## 3. Gauge manuscript correction — cross-repo (`~/projects/research/pmed-modern`), needs go-ahead

Every coverage number, the fold-split explanation, and the weak-identification
argument in `01-gauge-pmed/manuscript/gauge-pmed.qmd` come from the pre-fix
grids. Not a doc edit; a rerun.

- **Handoff first** (this repo can write it; the manuscript edit is the user's):
  the list of claims to change, with the post-fix numbers already measured
  (PR #34 / #35 FINDINGS) and what still needs the manuscript's own grid.
- **Rerun** (corrected 2026-08-22 night: Table 1 came from the hopper boot grid
  `inst/sim/hopper/run_gauge_boot_grid.R` running the *installed* pre-fix probmed,
  not from the local standalone): install the fixed dev package on hopper
  (`~/Rlib/4.4-postfix`, known-answer probe against the old lib), resubmit the same
  64-task array from `~/gauge_boot_postfix/` (pilot `--array=1,33` first, keep it,
  submit the complement), then `submit_collate.sh`.
- **Expected shape of the correction:** the Simulation section's "Wald
  under-covers, percentile conservative" story becomes "Wald calibrated,
  percentile no better"; the Application section's near-null / Fieller argument
  stands; the weak-ID diagnostic is demoted to a bootstrap-instability note or
  removed.

## 4. dev -> main release — needs a decision (package ON HOLD)

The online docs still describe the buggy estimator. A release is the only path
that updates them (pkgdown deploys on push to main). Version bump off `.9000`,
NEWS already written, `cran-comments.md` resync, `R CMD check --as-cran` against
CRAN medfit. ~1 h mechanical once the hold is lifted.

## Order

1 -> 2 in this session (no decisions needed). 3's handoff can be written now;
its rerun and the edit wait for the user. 4 waits for the hold.
