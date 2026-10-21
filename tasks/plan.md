# Implementation Plan: probmed single-mediator inference repair + open items

> **Superseded 2026-10-15:** This historical plan is not the active execution index. Its repair-plan/spec links point to files absent from this checkout, and its old checklist contains stale statuses. Use `tasks/plan-mc-bootstrap.md` and `tasks/todo.md` for active Phase 7 work; `.STATUS` owns current next actions. This file is retained as historical context.

> **Status 2026-10-14:** Phases 0–5 are complete (PRs #41–#45 merged). Phase 6 (T20–T22, `tasks/todo.md`) replaces the original "scale n_sim with N" follow-up with common random numbers, after a code review showed the simulation draws are fresh in every bootstrap draw.

Created 2026-10-13 (dev @ `4901862`). Task checklist: [`tasks/todo.md`](todo.md).

## Overview

This is the execution index for everything still open in probmed. It does **not**
restate the detailed requirements. Those live in the Codex repair plan and spec,
which stay authoritative for the technical content:

- `planning/repair/PLAN-probmed-repair-2026-10-09.md` (on branch
  `codex/probmed-repair-2026-10-09`): the ordered work, acceptance gates and
  author decisions of 2026-10-09
- `planning/repair/SPEC-probmed-repair-2026-10-09.md` (same branch):
  tolerances, the calibration matrix and the decision rule

This file adds three things: the steps needed to land the existing batches 1–2,
the gaps the 2026-10-13 review found, and the unrelated items that are still open.

### Current state (verified 2026-10-13)

| Item | State |
|---|---|
| Repair batches 1–2 | Branch `codex/probmed-repair-2026-10-09` (`3ac3de3`, `a1e3737`), base `5c6d837`. **Local only, not pushed.** Its worktree directory is gone. Full suite on an exported copy of the branch: **621 passed / 0 failed / 0 errors / 0 skips** |
| Overlap with dev since the base | `.Rbuildignore`, `AGENTS.md`, `scripts/check-local.sh` only, so the rebase should be clean |
| Defects still on dev | Non-Gaussian families are always simulated as Bernoulli (`R/compute-core.R:119,131`). The core ignores covariates. Bootstrap calls to `.pmed_core_simple` don't forward `n_sim` (`R/compute-bootstrap.R:42,135`). The parametric bootstrap holds residual SDs fixed |
| #38 gauge Table 1 rerun | Launched on hopper 2026-08-22. No post-fix gauge CSV on dev (`gauge_boot_coverage_nsim2000.csv` dates from 2026-07-16). Nothing posted to the issue |
| `.STATUS` `next:` item (1) | Stale: says "bump → 0.3.0", but 0.3.0 is already released |
| `AGENTS.md` (`4901862`) | Committed without the user's decision. Contains an absolute `/Users/davood/.codex/...` path, and the repo is **public** |

## Architecture Decisions

- **Land batches 1–2 before starting batch 3.** They are validated (621/0) and
  independent of the rest; holding them back only adds rebase risk.
- **Batch 3+ runs in sequence, not in parallel.** The family guard, `n_sim`/seed
  forwarding, standardization and the response bootstrap all touch the same
  `.pmed_core_simple` call sites (`compute-core.R`, `compute-bootstrap.R`,
  `compute-parallel.R`).
- **Fail early, then add capability.** The family guard comes before
  standardization, so unsupported models stop producing silently wrong numbers
  right away.
- **The response bootstrap is opt-in.** The default bootstrap stays unchanged
  until calibration passes (author decision, 2026-10-09).
- **Bug repairs may merge even if calibration fails.** Only calibrated-inference
  claims are restricted (author decision, 2026-10-09).
- **The submission/release hold stays.** No version bump and no dev→main.
- **All code goes on `feature/*` branches.** Each task names its target branch
  as a requirement. Creating or switching a branch still needs the user's
  explicit go-ahead each time.

## Dependency graph

```
T1 back up branch ─► T2 clean + rebase ─► T3 NA-rate/runtime measurement ─► T4 R CMD check ─► T5 PR (batches 1–2)
                                                │
                                                └─(if NA rate high)─► T3b basin-start fix ─► T4
T5 merged ─► T6 family guard ─► T7 n_sim/seed forwarding ─► T8 covariance validation
                                                        └─► T9 standardization (Gaussian) ─► T10 standardization (binary)
T10 + T11 separation detector + T8 ─► T12 response bootstrap: simulate/refit ─► T13 diagnostics ─► T14 docs
T14 ─► T15 cost benchmark ─► (user approval) ─► T16 calibration pilot on hopper ─► T17 full calibration
T18 #38 rerun collation ─ independent of all of the above
```

## Task List

Full acceptance criteria and verification steps for each task are in `tasks/todo.md`.

### Phase 0: Secure the work (needs the user's yes)
- [ ] T1: Push `codex/probmed-repair-2026-10-09` to origin as a backup

### Phase 1: Land batches 1–2 (branch `feature/mbco-repair`)
- [ ] T2: Rebase on dev, rename to `feature/*`, strip the manuscript planning files and the `.STATUS` line
- [ ] T3: Measure the MBCO NA rate and runtime (decision point)
- [ ] T3b: *(only if T3 fails its threshold)* More robust wrong-sign-basin starts
- [ ] T4: `R CMD check --as-cran` with the roxygen toolchain gate
- [ ] T5: Open a PR to dev with evidence; author reviews the ML-variance point-estimate change

### Checkpoint A: batches 1–2 merged
- [ ] dev full suite green; CI green; the author has accepted the point-estimate convention change

### Phase 2: Support boundaries and stochastic controls (branch `feature/repair-batch3`)
- [ ] T6: Reject unsupported family/link/structure before any random draws
- [ ] T7: Validate and forward `n_sim` and `seed` on every stochastic path (single + parallel)
- [ ] T8: Covariance-extraction validation against an independent stacked reference

### Checkpoint B
- [ ] Unsupported models fail early with informative errors; `n_sim`/`seed` documented and tested; covariance disposition recorded

### Phase 3: Covariate standardization (branch `feature/repair-standardize`)
- [ ] T9: Observed-sample standardization, Gaussian path (cancellation preserved)
- [ ] T10: Standardization, binary path, validated against independent g-computation
- [ ] T11: Binary separation detector (original fit and per-draw)

### Checkpoint C
- [ ] Adjusted binary example matches the g-computation reference; separation detector fires on a planted separable fixture

### Phase 4: Opt-in response bootstrap (branch `feature/repair-response-boot`)
- [ ] T12: Simulate M→Y, refit both equations on the same rows, recompute standardized P_med
- [ ] T13: Completion and failure diagnostics (separation, non-convergence) on the result
- [ ] T14: Docs: roxygen, README.qmd→md, NEWS, vignette

### Checkpoint D
- [ ] Full suite plus `R CMD check` green; opt-in path documented; default unchanged

### Phase 5: Verification and calibration (needs user approval for compute)
- [ ] T15: Cost benchmark on the spec's calibration matrix
- [ ] T16: Hopper pilot via `sbatch --array=1-1`, then the full array
- [ ] T17: Collate PASS/FAIL/INCONCLUSIVE per cell; restrict claims where needed

### Independent track
- [ ] T18: #38: collate the post-fix gauge rerun from hopper, diff against Table 1, draft the issue comment (posting needs a yes)
- [ ] T19: Housekeeping: `.STATUS` stale item (1) and repair status on the next real `.STATUS` edit; the `AGENTS.md` decision

### Checkpoint: Complete
- [ ] All six acceptance gates in the Codex PLAN are met or have an explicit disposition
- [ ] The author has reviewed every change to inference, numerical results or the bootstrap target before manuscript claims are revised

## Risks and Mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| Strict "any grid point failed → interval NA" makes MBCO intervals missing too often | High | T3 measures the rate before the PR; T3b if it is above threshold |
| Grid size (~200 P_med points, ~960 IE points × 2 sign cases × starts) makes calibration expensive | Med | T15 cost benchmark before any hopper run; consider a coarser IE scan plus local refinement |
| Codex regenerated `man/` with roxygen2 8.1.0, and 8.x breaks the docs for S7 methods (memory `roxygen-version-toolchain`) | Med | T4 gate: no codoc or `\usage` warnings; otherwise re-run `document()` under 7.3.3 |
| The ML-variance (SSE/n) point-estimate change shifts `mbco` numbers in the manuscript | Med | Author review in T5; manuscript changes stay a separate, pmed-side task |
| medfit drift: installed 0.4.0, source 0.6.0, native SEM planned for 0.7.0 | Med | T8 tests both available versions; don't assume native SEM (Codex PLAN, "medfit dependency update") |
| The bounded scan misses disconnected components narrower than its spacing | Low | Already documented; `status = "ok"` is not a connectedness proof |
| Parallel MBCO still uses the old inversion (failed fit counted as rejection) | Med | Open question Q3 |

## Resolved defaults (2026-10-14, `/goal` "finish all the plans and commit and push and sync")

1. **T1:** the goal authorizes pushing. The backup is the *cleaned* `feature/mbco-repair`,
   not the raw Codex branch: the repo is public, and the raw branch carries pmed-owned
   planning docs and a handoff JSON.
2. **T2:** drop `planning/repair/*` and the handoff JSON from the package branch. pmed
   already tracks its own copies (`docs/PLAN-…`, `docs/SPEC-…`,
   `docs/reports/probmed-repair-handoff-2026-10-09.json`), so nothing is written to pmed.
3. **Parallel MBCO:** out of scope for this repair, and stated as a limitation in NEWS.
4. **T3 threshold:** MBCO NA rate ≤ 1% on well-conditioned cells, ≤ 5% near the null.
5. **`AGENTS.md`:** keep `4901862` as committed (matches medfit `b2c3c4a`).

Still needs the user's explicit yes: merging any PR, the full hopper calibration array
(after the pilot and cost projection), and posting to #38.
