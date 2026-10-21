# Task list: probmed CI tiers

> **Superseded 2026-10-20:** Retained only as the checklist for the historical
> [hosted-tier proposal](plan-ci-tiers.md). The adopted policy and safe switch order
> are in [the local-CI transition plan](PLAN-local-ci-transition-2026-10-20.md).
> Unchecked tasks below are not an active implementation queue or evidence of
> completed work. Preserve the scientific merge hold and exact-head local gate.

Track this work independently from the unfinished general plan in `tasks/todo.md`.
Scope: R package repositories only; these tasks target probmed. Any other R
package repository needs its own branch-policy and CI audit before adopting this tiering.

## Phase 1: Live constraints

### T1: Audit branch rules and current CI cost

**Description:** Verify `dev`/`main` branch protection, merge strategy, required status contexts, active branch names, and recent run durations.
**Acceptance criteria:**

- [ ] Record whether `dev` and `main` are PR-gated and list exact required check names.
- [ ] Record queue time, execution time, and runner minutes for recent R-CMD-check, coverage, and pkgdown runs.
- [ ] Confirm whether `master` remains active and whether pushes bypass the PR path.
**Verification:** Inspect live repository settings and workflow run details; do not rely solely on stale `CLAUDE.md` notes.
**Dependencies:** None · **Scope:** S

## Checkpoint A: Policy confirmed

- [ ] If `dev` is not PR-gated, retain full remote checks on `dev` pushes or revise branch enforcement before continuing.

## Phase 2: Gate contract

### T2: Specify fast and full gate contents

**Description:** Define what runs on feature → `dev` versus `dev` → `main`, including tests skipped by `NOT_CRAN`, coverage, platform matrix, package docs/vignettes, pkgdown, and stable status reporting.
**Acceptance criteria:**

- [ ] Fast feature PR gate is a single Ubuntu/current-R remote job and catches deterministic package failures.
- [ ] Full main PR gate retains the existing cross-platform R-CMD-check, coverage, and pkgdown validation.
- [ ] Every required PR gate reports a stable status and cannot silently remain pending because a job was skipped.
**Verification:** Review trigger/condition table against live branch rules and package test conventions.
**Dependencies:** T1 · **Scope:** S

## Phase 3: Local preflight clarity

### T3: Align local command guidance with actual remote CI

**Description:** Correct `scripts/check-local.sh` and contributor docs so local preflight is not described as equivalent to current PR CI.
**Acceptance criteria:**

- [ ] Explain the local checks, deliberate skips, logs, and invocation.
- [ ] State that local hooks/preflight are advisory and GitHub required checks enforce merges.
- [ ] Keep automatic hooks out of scope unless separately authorized.
**Verification:** Run documented commands or inspect their current behavior; review all edited guidance for stale equivalence claims.
**Dependencies:** T2 · **Scope:** S

## Checkpoint B: Implementation ready

- [ ] Gate contract and docs are agreed before workflow behavior changes.

## Phase 4: Workflow implementation

### T4: Implement branch-specific remote CI tiers

**Description:** Configure fast feature → `dev` CI and full `dev` → `main` CI, removing only verified duplicate push runs.
**Acceptance criteria:**

- [ ] Feature PR to `dev` runs and reports the fast required remote gate.
- [ ] PR to `main` runs full R-CMD-check matrix, coverage, and pkgdown validation.
- [ ] Push triggers are removed only where PR-only branch updates are confirmed; pkgdown deploy on `main` remains intact.
**Verification:** Validate YAML/conditions and required check reporting; use actionlint/repository workflow checks if available.
**Dependencies:** T1, T2, T3 · **Scope:** M

## Phase 5: End-to-end verification

### T5: Verify runtime, required checks, and branch integration

**Description:** Exercise both PR routes and compare runtime/cost with the baseline.
**Authorization gate:** Obtain explicit user approval before opening the validation PR; obtain separate explicit approval before merging.
**Acceptance criteria:**

- [ ] Feature → `dev` fast gate passes with no missing or pending required status.
- [ ] `dev` → `main` full gates pass and remain enforceable.
- [ ] Report before/after wall time and runner-minute change; reconcile project status/docs.
**Verification:** Inspect live check results and branch settings after each workflow run.
**Dependencies:** T4 · **Scope:** S–M

## Completion gate

- [ ] Both PR paths have the intended remote gates and stable required statuses.
- [ ] Local preflight is clearly advisory.
- [ ] No required check is skipped/pending, and full platform/coverage validation still precedes main.
- [ ] Current project status and working tree are reconciled.
