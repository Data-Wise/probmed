# Boundary point-risk: task list

Plan: [`plan-boundary-point-risk.md`](plan-boundary-point-risk.md). Created 2026-10-21.
Source of the failing evidence: PR #55 branch `codex/dev-main-review`,
`docs/reports/VALIDATE-wasserstein-boundary-selection-2026-10-20.md` and
`docs/reports/data/wasserstein-selection-2026-10-20/replicates.csv`.
Research checks: `Rscript scripts/test-wasserstein-boundary-selection.R`,
`Rscript scripts/audit-wasserstein-boundary-selection.R`.

## Phase A: understand the failure

- [x] **T1 Risk/bias derivation** (S, 1 doc) — done 2026-10-21, `docs/specs/DERIVE-boundary-point-risk-2026-10-21.md`
  - Accept: derives bias and risk of an axis-snapping estimator as a function of the
    snap probability and the distance of the true score from the axis, for the regular,
    zero-direct, zero-indirect and weak-direct cells; states what any threshold rule
    must trade off.
  - Verify: every claim is derived or cited to a source actually read; no claim about a
    paper not opened.
  - Files: `docs/specs/DERIVE-boundary-point-risk-<date>.md`
  - Deps: none

- [ ] **T2 Check against archived replicates** (S, 1 script + output) — RUN 2026-10-21, overall FAIL on the SD gate (bias 20/20, snap 10/10, SD 17/20; both controls rejected as required). Deviations from the criteria below: covers the 5 binary cells only (cells 1, 3 are Gaussian), and uses a seeded Gaussian-experiment prediction, so it does use RNG (seed 91, not a study seed). Follow-up (`scripts/check-boundary-risk-studentization.R`): shuffling v_hat reproduces the cell 5 and 7 select995 misses (studentization explains 2 of 3); cell 4 select95 miss traced to d_hat tail shape/scale (`scripts/check-boundary-risk-cell4.R`; not sampling noise, P = .003; v_hat variability irrelevant); pair bootstrap passes 10/10 but cannot fail. Left unchecked; author decides at Checkpoint A.
  - Accept: the T1 formulas, evaluated with the observed snap rates (22.95%, 52.70%,
    94.95%, 99.65%) and cell truths, reproduce the reported bias/SD within stated
    tolerance in all 7 cells; includes a positive control (a cell where the formula must
    fail if wrongly specified).
  - Verify: script run, output quoted; no new simulation or RNG.
  - Files: `scripts/check-boundary-risk-derivation.R`
  - Deps: T1

### Review A: Codex adversarial review (after T2)
- [x] Run `/codex:adversarial-review` on `feature/boundary-point-risk` (verdict: needs-attention, 2 high, 1 medium)
- [x] Recorded verbatim with triage in `docs/specs/REVIEW-codex-boundary-point-risk-A-2026-10-21.md`
- [ ] No unresolved high-severity finding: finding 1 fixed (D3 corrected); **finding 2 open**
      (T2 does not test the derivation) -> T2b below. Finding 3 fixed by rewording.

- [ ] **T2b Test the derivation itself** (S, 1 script + output) — RUN 2026-10-21: A1-A3 pass (slope .470 (corrected from .459, wrong index), coefficient .964, control 1.913), B1-B3 10/10, **C1 8/10 (FAIL)**: cell 4 select95 and cell 5 select995, the known rare-snap cells. Tolerances fixed before the run, not changed. Left unchecked pending Checkpoint A.
  - Accept: (a) a rate check: analytic or large-draw bias of the raw plug-in at an exact
    zero contrast across a grid of s (equivalently n), log-log slope near 1/2 with a stated
    tolerance, plus a planted wrong-rate control; (b) D2's snap formula computed in closed
    form (independence product and, separately, with the correlation) compared with
    simulation per cell, stated tolerance; (c) a squared-risk decomposition
    (bias^2 + variance, with conditional snap terms) reproducing the observed RMSE per cell.
  - Verify: script run, output quoted, tolerances fixed before the run; failures reported.
  - Files: `scripts/check-boundary-risk-rate.R`, `docs/specs/data/`
  - Deps: T2, Review A

### Checkpoint A (author) — ACCEPTED 2026-10-21
- [x] Derivation matches the observed bias pattern in the binary cells; SD misses have
      supported (not proven) explanations. T2 verdict FAIL 17/20 stands.
- [x] Review A triaged (finding 2 addressed by T2b)
- [x] Author approves continuing to Phase B, with these limitations stated:
      (1) T2b C1 fails 8/10 in the rare-snap cells (cell 4 select95, cell 5 select995);
      (2) D3 (no threshold improves the n^(-1/4) rate) was an argument at acceptance; now a local-minimax two-point lemma D3' pending independent audit;
      (3) the cell 4 mechanism is a hypothesis, not identified.

## Phase G: gate design (added 2026-10-21; blocks T4/T5)

Source: `docs/specs/GRILL-boundary-gates-2026-10-21.md`, `docs/specs/REVIEW-codex-boundary-point-risk-G-2026-10-21.md`.

- [ ] **G1 Author writes per-cell gates** (author; XS) — HELD on the author
  - Accept: per-cell risk/bias gates stated on substantive grounds, separate tolerances at n = 500 and n = 1000 (Q5 reversal), no reference to any candidate's archived risk.
  - Deps: none
- [x] **G2 Floor and ladder derivation** (S, 1 doc + script) — DONE 2026-10-21 at n = 500 only: floor on max RMSE about .09 (direct ladder) and .07 (indirect); `docs/specs/DERIVE-boundary-floor-2026-10-21.md`. n = 1000 and ladder-to-ladder pairs still open.
  - Accept: two-point/minimax floor per ladder point from the DGP laws alone; ladder positions set by data-law TV along the DGP path (not SE units of d2); zero-indirect side (cell 6 analog) included.
  - Verify: script generalizes `scripts/check-gate-feasibility-tv.R`; controls (shift 0, large shift) reproduced; output saved under `docs/specs/data/`.
  - Files: `docs/specs/DERIVE-boundary-floor-<date>.md`, `scripts/check-boundary-floor-ladder.R`
  - Deps: none (can run before G1; available on request)
- [ ] **G3 Feasibility pre-check of the author's gates** (S) — TOOL READY 2026-10-21: `Rscript scripts/check-gate-tolerance.R <n> <tol_direct> [<tol_indirect>]` (controls pass; verdict INFEASIBLE / NEAR-FLOOR / ROOM against the G2 floor). Left unchecked: needs the author's tolerances (G1). The 1.5x NEAR-FLOOR cut is a reporting convention, not derived.
  - Accept: each proposed gate is tested against the G2 floor; any infeasible gate is reported as a finding, not retuned.
  - Deps: G1, G2
- [ ] **G4 One complete gate set** (S, plan amendment) — SKELETON DRAFTED 2026-10-21 (`docs/specs/SPEC-boundary-gate-set-DRAFT-2026-10-21.md`); all AUTHOR cells blank, not registered
  - Accept: ladder, tolerances per n, boundary flag definition and miss-rate gate, between-ladder-point risk control, and which of bias/SD, SE/SD, coverage are kept; Q1-Q3, Q6 promoted from PROVISIONAL or reversed in the ledger.
  - Deps: G3
- [ ] **Review G2 (Codex)** on the complete gate set, before any freeze (this is the Review B focus for gates; run once, record verbatim, triage)
  - Deps: G4

### Checkpoint G (author)
- [ ] Gate set registered in a commit before any fresh data; feasibility pre-check passed or the negative-result endpoint chosen

## Housekeeping (not blocking)
- [ ] `.gitignore` entry for private teaching tarballs (`wasserstein-inference-private-*.tar.gz`, ~73 MB) in the main checkout
- [ ] Check the external `AGENTS.md` reference under `~/.codex/reference/` exists
- [ ] F7 teaching files (pmed-modern, uncommitted by instruction) need an Obsidian vault copy to survive re-sync; commit only on the user's go
- [ ] Reconcile T2/T2b boxes after the author accepts the recorded FAILs as results (currently left unchecked on purpose)

## Phase B: candidate and tests

- [x] **T3 Candidate specification** (M, 1-2 docs) — DONE 2026-10-21, `docs/specs/SPEC-boundary-candidate-T3-2026-10-21.md`. Result: smooth family C(kappa, tau) specified and predicted to fail the bias gate in cells 5, 6 at every grid point; **gate feasibility check shows no [0,1]-valued estimator can pass bias/SD < .10 in cells 5 and 7 together at n = 500 (TV .141 < needed .168-.179)**. T4/T5 paused pending author decision on the gate set; nothing frozen.
  - Accept: one smooth constrained estimator family (no hard threshold), its
    uncertainty estimate with an explicit justification, and written pass/fail
    expectations for regular and weak-signal cells; says which of the derived
    trade-offs it accepts.
  - Verify: a reviewer can predict, from the doc alone, in which cells it should fail.
  - Files: `docs/specs/PLAN-wasserstein-smooth-boundary-<date>.md`
  - Deps: T2

- [ ] **T4 Independent numerical tests** (M, 2-3 files)
  - Accept: scalar vs vectorized agreement, scale invariance, axis/symmetry cases,
    independent high-draw reference for the SE proxy within 10%, and a planted-defect
    check that the tests fail on a wrong implementation.
  - Verify: `Rscript scripts/test-wasserstein-smooth-boundary.R` pass counts quoted;
    planted defect caught.
  - Files: `R/` or `scripts/` candidate code, `scripts/test-...R`
  - Deps: T3, Checkpoint G (a candidate cannot be tested against an infeasible gate)

### HOLD 2026-10-21: T4 and T5 are held until Checkpoint G (author's per-cell gates, G1-G4).
### Codex review G of the infeasibility argument is done and triaged.

### Review B: Codex adversarial review (after T4, BEFORE T5)
- [ ] (Deps: Checkpoint G, T4) Run `/codex:adversarial-review` with focus: is the design the right one, which
      assumptions can fail, what could be tuned after the fact (gate, seed, tolerance,
      estimator), can the tests fail
- [ ] Record verbatim in `docs/specs/REVIEW-codex-boundary-point-risk-B-<date>.md`, triage
      every finding; accepted findings are applied by us, then re-review only if the design changed
- [ ] No unresolved high-severity finding before T5

- [ ] **T5 Freeze commit** (XS) — only after Review B is triaged; records the triage
  - Accept: design, estimator, tolerances, gates (unchanged), seeds (disjoint from all
    prior studies, audited), and source SHA-256 committed before any fresh dataset.
  - Verify: `git log` shows freeze commit; seed-disjointness audit output quoted.
  - Deps: T4

### Checkpoint B (author, go/no-go)
- [ ] Tests pass, nothing about the gates changed, no fresh data generated yet
- [ ] Review B triaged and the frozen commit matches the reviewed design
- [ ] Author explicitly approves the fresh-seed run (T6)

## Phase C: validation (not started without approval)

- [ ] **T6 Fresh-seed validation** (M)
  - Accept: unchanged gates on all 7 cells, zero failures, MCSE, selection rates, raw
    rows and provenance archived.
  - Verify: independent audit script reproduces every summary row.
  - Deps: Checkpoint B

### Review C: Codex adversarial review (after T6, before the T7 verdict)
- [ ] Run `/codex:adversarial-review` with focus: retuning or post-hoc gate reading, seed
      disjointness, provenance/hashes, whether the audit reproduces every summary row,
      whether the verdict overclaims
- [ ] Record verbatim in `docs/specs/REVIEW-codex-boundary-point-risk-C-<date>.md`, triage
      every finding
- [ ] No unresolved high-severity finding before T7

- [ ] **T7 Decision report and doc sync** (S) — includes the Review C triage
  - Accept: pass/fail verdict per cell; PRs #54/#55 stay blocked unless every unchanged
    gate passes; `.STATUS` updated once, without restating the report.
  - Deps: T6

## Not in scope

I7, I6 reruns, production adoption, merge, required-status switch, release, submission.
