# Implementation Plan: boundary point-estimation risk / bias trade-off

Created 2026-10-21 (dev @ `eb6516f`). Task checklist: [`todo-boundary-point-risk.md`](todo-boundary-point-risk.md).
Does not touch the older `tasks/plan.md` / `tasks/todo.md` (superseded Phase 0-6 work, left as is).

## Status (2026-10-21)

Phase A and T3 done; Reviews A and B1 triaged; grill ledger Q1-Q6 recorded (Q5 reversed, Q1-Q3 and Q6 provisional). T4/T5 held pending the author's per-cell gates. See `docs/specs/` and `.STATUS`.

## Overview

`.STATUS` next research action: derive the boundary point-estimation risk/bias
trade-off before freezing a new candidate (PRs #54/#55 stay draft until the
unchanged complete diagnostics pass). Evidence so far
(`docs/reports/VALIDATE-wasserstein-boundary-selection-2026-10-20.md`, PR #55 branch
`codex/dev-main-review`): Fieller intervals pass coverage in every binary cell, but
neither hard-selection point estimator (lambda = qnorm(.975), qnorm(.9975)) passes any
binary cell's bias/SD, SE/SD and coverage gates together. Hard selection snaps to an
axis in 22.95% (regular n=500) and 52.70% (strong) of fits, which induces bias where the
true score is interior; strong selection also leaves the point outside its own interval
in up to 30% of fits.

This plan produces a **derivation, a frozen candidate design, and independent numerical
tests**. It stops before any fresh-seed simulation and before any production change;
those need the author's approval.

## Hard constraints (from `.STATUS`, PR #54/#55 body)

- No gate relaxation, no threshold escalation to chase exact-null bias, no change to any
  frozen diagnostic, I7 closed, no merge / required-status switch / release / submission.
- Preserve all regular and weak-signal cells; seed ranges disjoint from every earlier
  study (earlier: 20561020 + cell*10000 + replicate and its predecessors).
- Design frozen (committed) before any fresh dataset is generated.
- Research code lives on `codex/dev-main-review` (PR #55). New work goes on a new
  feature branch cut from it; **creating that branch/worktree needs the user's explicit
  instruction** (no-unrequested-branch-switch rule).

## Dependency graph

```
T1 risk/bias derivation (math, per cell)
  ├── T2 numerical check of derivation against existing replicates.csv (no new sims)
  │     └── T3 candidate family + frozen design (smooth constrained estimator + uncertainty)
  │           ├── T4 independent numerical tests (scalar/vector, scale, axis, precision)
  │           └── T5 pre-registered freeze commit (design + source hash + seeds)
  │                 └── T6 fresh-seed validation  [AUTHOR GATE]
  │                       └── T7 decision report + .STATUS reconciliation
  Codex reviews: R-A after T2 | R-B after T4, before T5 | R-C after T6, before T7
```

Slices are vertical: each task yields a reviewable artifact (derivation note, check script
with output, frozen design, tests) rather than a layer.

## Phases

### Phase A: Understand the failure (read/derive only)
- T1 Derive risk/bias trade-off for axis-snapping estimators
- T2 Check the derivation against the archived replicates (reuses data, no new RNG)

### Review A: Codex adversarial review (after T2, before Checkpoint A)
- See "Codex review gates" below. Focus: does the derivation survive scrutiny, and do the
  T2 checks (and their follow-ups) actually test it or only restate the archive.

### Checkpoint A (author review)
- Derivation predicts the observed bias pattern across the binary cells and every SD miss
  has an identified mechanism (met 2026-10-21; T2 verdict stays FAIL 17/20 as recorded).
- Review A triaged; no unresolved high-severity finding.

### Phase B: Candidate and tests
- T3 Specify smooth constrained candidate (no hard threshold) and its SE. Constraints from
  Phase A: SD/SE acceptance only by fresh-DGP simulation (Gaussian-experiment SD is
  unreliable in rare-snap cells); carry the joint (d_hat, v_hat) dependence; check the
  literature on unbiased estimation of sqrt|delta| before citing it.
- T4 Independent numerical tests for the candidate

### Review B: Codex adversarial review (after T4, BEFORE the T5 freeze commit)
- Must precede the freeze: after T5 nothing may change. Focus: is the design the right
  one, which assumptions can fail, which frozen element (gate, seed, tolerance, estimator)
  could be tuned after the fact, are the tests able to fail.

### Phase B (cont.)
- T5 Frozen design commit, only after Review B is triaged; the freeze commit records the
  Review B triage.

### Checkpoint B (author review, go/no-go for fresh seeds)
- Design, thresholds, tolerances, seeds, estimators frozen; tests pass; no gate changed.
- Review B triaged and the frozen commit matches the reviewed design.

### Phase C: Validation (only after approval)
- T6 Fresh-seed run, unchanged gates, all rows archived

### Review C: Codex adversarial review (after T6, before the T7 verdict)
- Focus: any retuning or post-hoc gate reading, seed disjointness, provenance and hashes,
  whether the independent audit reproduces every summary row, and whether the verdict
  claims more than the evidence supports.

### Phase C (cont.)
- T7 Report, `.STATUS`/`CLAUDE.md` sync; PRs stay blocked unless every unchanged gate passes.
  The report includes the Review C triage.

## Codex review gates

One adversarial review per phase (A, B, C), run with `/codex:adversarial-review` against
the feature worktree (`feature/boundary-point-risk`), with the phase focus text above.
Rules:

- **Review-only.** Codex output is recorded verbatim in
  `docs/specs/REVIEW-codex-boundary-point-risk-<phase>-<date>.md`; it applies no fixes.
- **Triage, do not auto-apply.** Each finding gets accept / reject / defer with a reason,
  in a table in the same file. Accepted findings become tasks or edits made by us, not by
  the reviewer. A rejected finding needs a stated reason, not silence.
- **Blocking.** An unresolved high-severity finding blocks the checkpoint that follows.
  A "no material findings" verdict is recorded but is not validation: an earlier review of
  an almost-empty working tree returned "approve", which says nothing about the science.
- **Scope honesty.** Codex sees the working tree or branch diff. The review states which
  files and data it could not see (for example the archived CSV on the PR #55 branch).
- **Order matters for Phase B.** Review B runs before the T5 freeze commit, not after.
- **No gate relief.** A review cannot relax a frozen gate or authorize the T6 run; only the
  author does that at Checkpoint B.
- **Cost and cadence.** Three reviews total; no review per task. Re-run a phase review only
  if its accepted findings change the design.

## Risks

| Risk | Impact | Mitigation |
|------|--------|------------|
| Derivation explains bias but implies no estimator passes all gates | High | Say so in T1/T2; a negative result is a valid outcome and is reported |
| Candidate tuned after seeing results | High | T5 freeze commit precedes any fresh data; no retuning |
| Smooth estimator passes exact-null but fails weak-signal (earlier pattern) | High | T3 criteria include weak-direct and regular cells up front |
| SE from Gaussian-proxy again miscalibrated (SE/SD 0.49-2.8) | High | T3 states the SE justification; T4 checks proxy precision independently |
| Work lands on the wrong branch / collides with PR #55 session | Med | Branch is `feature/boundary-point-risk`, cut from `dev` (author decision 2026-10-21); archived data read from the PR #55 checkout |
| Codex review treated as validation | Med | Review gates are triage inputs only; "approve" is recorded, not credited; no gate relief |
| Review B after the freeze makes its findings unusable | High | Review B is sequenced before T5 |

## Open questions for the author

1. (Resolved 2026-10-21) Branch from `dev`.
2. Is a negative result (no admissible estimator) an acceptable endpoint for this plan?
3. Budget for T6 (Hopper vs local; earlier study took 318 s locally for 44,000 rows).
