# Boundary point-risk: decision report (2026-10-21)

Branch `feature/boundary-point-risk` (from `dev`), plan `tasks/plan-boundary-point-risk.md`. Endpoint chosen by the author on 2026-10-21 ("option 1"): the
negative result. No production code changed; no frozen diagnostic changed; PRs #54/#55 stay held under their original gates (Q4).

## Verdict

No point estimator tested, and no estimator that the frozen gate could ever accept, resolves the boundary point-estimation problem for the binary-outcome
Wasserstein proportion mediated. The supported deliverable remains the Fieller interval, which passes coverage in every archived binary cell. T4 (tests), T5 (freeze),
T6 (fresh-seed validation), Review B/C and the T7 gate run are NOT carried out: the pre-committed pilot rule (grill Q8) stops them.

## Findings (each with its evidence file)

| # | Finding | Status | Evidence |
|---|---------|--------|----------|
| 1 | Plug-in bias at an exact zero contrast is of order n^(-1/4), coefficient E\|Z\|^(1/2) = .822; the "leading-order" label is a ratio-of-means heuristic | derived, verified (T2b rate slope .470) | `DERIVE-boundary-point-risk-2026-10-21.md` D1; `data/` |
| 2 | Snap/threshold estimators move the bias; no threshold improves the rate (local-minimax two-point lemma D3') | derived, pending independent audit | DERIVE D2, D3' |
| 3 | The frozen bias/SD < .10 gate is infeasible for any [0,1]-valued estimator in cells 5 and 7 at n = 500 (TV .1404 vs needed .179); empirical-pass probability <= 8.3e-5 | audited | `SPEC-boundary-candidate-T3-2026-10-21.md` section 2; `check-gate-feasibility-tv.R` |
| 4 | Two-point floor on max RMSE along the weak-signal ladders: about .09 (direct) and .07 (indirect) at n = 500; .08 and .06 at n = 1000; the boundary pair binds | computed, controls pass; multipoint bound not computed | `DERIVE-boundary-floor-2026-10-21.md` |
| 5 | A gate set v1.1 was drafted and reviewed (Codex G2; 5 findings applied): RMSE over 9-point ladders, unconditional 95% coverage, SE/SD on regular cells, boundary flag | proposed, never registered | `SPEC-boundary-gate-set-2026-10-21.md` |
| 6 | Pilot of five reference estimators (plug-in, C(kappa, .5), kappa in {0, 1, 2, 3.84}) against the author's tolerances (.15 / .12 / .125 / .09): all FAIL; best max RMSE .218 / .204 / .183 / .138, at the exact boundary | run, NEGATIVE | `data/pilot-g8-2026-10-21.txt`, `-by-point.txt` |

## What this does and does not show

- It shows that the old bias/SD gate cannot be met, that the risk floor and the trade-off make boundary bias and interior risk compete, and that the five reference estimators
  cannot meet the author's tolerances. The failures concentrate at the exact boundary (bias about -.21 direct, +.18 indirect for C(3.84)).
- It does not show that no estimator can meet the tolerances: the floor (.09, .07) is below them. A different family (for example a larger shrinkage or a flag-based hybrid)
  was not tested because continuing requires the author to name a new family before any further simulation.
- Known limits kept on record: T2 SD gate 17/20 and T2b C1 8/10 (rare-snap cells, studentization), D3' pending audit, multipoint floor uncomputed, flag tables are
  Gaussian approximations, and a small R = 100 smoke run preceded the full pilot (design committed first as `58029e7`).

## Reproduction

`scripts/check-boundary-risk-*.R` (T2/T2b), `scripts/check-gate-feasibility-tv.R`, `scripts/check-boundary-floor-ladder.R [n]`, `scripts/check-boundary-floor-pairs.R [n]`,
`scripts/check-gate-tolerance.R`, `scripts/decide-gate-rmse.R`, `scripts/audit-seed-ranges.R`, `scripts/tabulate-flag-misses.R`,
`scripts/pilot-reference-estimators.R <out_dir> 2000` (set `VECLIB_MAXIMUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1` first, or the forked workers stall on BLAS threads).
The audit and pilot scripts need the PR #55 checkout path (env `PROBMED_PR55_SCRIPTS` / `PROBMED_PR55_DIR`, or an argument); they stop with a message if it is missing.
Archived data (read only): `docs/reports/data/wasserstein-selection-2026-10-20/replicates.csv` on the PR #55 branch checkout.

## Held and open

PRs #54/#55 held under their original gates; any decision about an interval-only path is separate and explicit. Teaching copy `F7-*` in pmed-modern is uncommitted (needs a vault copy).
Reopening requires the author to name a new candidate family; Review G3 and Checkpoint G then apply to the revised gate set.
