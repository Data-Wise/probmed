# Implementation Plan: Tiered local and GitHub CI for probmed

> **Superseded 2026-10-20:** Historical proposal, not the adopted CI policy.
> The author selected commit-bound local validation with a locally posted required
> status instead of the hosted tiers proposed below. Use
> [the local-CI transition plan](PLAN-local-ci-transition-2026-10-20.md).
> Its protection switch remains gated on exact-head validation and the scientific
> merge hold. Do not implement this older proposal or treat its checklist as active.

## Scope

This policy applies only to R package repositories. This plan is specifically for
`Data-Wise/probmed`; applying the same pattern to another R package requires a
separate audit of that repository's branch model, tests, and CI. Do not apply it
to manuscript, book, or other non-package repositories.

## Overview

Reduce repeated long GitHub Actions runs while preserving a dependable remote merge gate at each integration boundary. Use a fast remote Linux gate for feature pull requests into `dev`; run the full platform matrix, coverage, and site validation on `dev` → `main`. Keep local preflight as earlier feedback, not as the merge authority. Do not change workflows until the current branch protections and required status checks are confirmed.

This is a standalone plan. The existing `tasks/plan.md` and `tasks/todo.md` contain unfinished, unrelated work and remain untouched. The task checklist is `tasks/todo-ci-tiers.md`.

## Evidence and current state

- `.github/workflows/R-CMD-check.yaml` runs five jobs (macOS release, Windows release, and Ubuntu devel/release/oldrel-1) on pull requests and pushes to `main`, `master`, and `dev`.
- `.github/workflows/test-coverage.yaml` runs on the same PR and push branches, so a PR merge into `dev` can repeat the full coverage job.
- `.github/workflows/pkgdown.yaml` validates PRs to `dev` without deployment; it deploys on non-PR events, including pushes to `main`/`master`.
- `scripts/check-local.sh` runs `devtools::test()` with `NOT_CRAN=false`, then builds/checks a tarball without tests or vignettes. Its comments describe it as matching PR CI, but the current PR workflows run a broader matrix and coverage; update that description as part of this work.
- PR #51 to `dev` ran R-CMD-check, test-coverage, and pkgdown remotely; all passed. The R-CMD-check matrix took about 23.5 minutes to finish. This is one observed PR, not a general runtime estimate.
- Current branch-protection settings could not be read through the GitHub integration. `CLAUDE.md` records a required Ubuntu release check on `main`, but that record must be verified against live settings before it is used to design required-check names.

## Architecture decisions

1. Preserve remote CI on feature → `dev`; `dev` is the shared integration branch, so local-only validation is not an adequate protection.
2. Use a single Linux/current-R PR gate on feature → `dev`; leave expensive simulation tests under the package's existing `NOT_CRAN` conventions and run the local fast preflight before push.
3. Run the complete five-job R-CMD-check matrix and full coverage on `dev` → `main`; retain pkgdown validation on PRs. Deployment is NOT only on the main push: `pkgdown.yaml` also triggers on `release` and `workflow_dispatch`, and its deploy step runs for every non-PR event, so a manual dispatch from any ref can publish. Phase 1 must audit these routes and restrict deploy to the main ref (or document and gate each route).
4. Remove duplicate R-CMD-check and coverage runs on pushes to `dev`/`main` only after confirming both branches accept changes through PRs and the required PR statuses block merging. Keep a post-merge smoke check only if the audit finds a concrete need.
5. Avoid path filters in the first implementation. If docs-only optimization is later warranted, keep an always-reported aggregate required status so a skipped workflow cannot leave a required check pending.
6. Keep `scripts/check-local.sh` advisory. A local hook may be absent or bypassed and does not replace GitHub's remote check.

## Dependency graph

```text
T1 live policy and duration audit
 └─> T2 define branch-specific gate contract and stable required status names
      └─> T3 reconcile local preflight documentation
           └─> T4 implement workflow tiers and remove verified duplicate runs
                └─> T5 exercise feature→dev and dev→main paths; review and merge
```

## Phases and checkpoints

### Phase 1: Confirm live constraints

- Audit branch protections/rulesets for `dev` and `main`, merge queue settings, required status contexts, and whether direct pushes are possible.
- Collect elapsed time and runner-minute data from several recent PRs for R-CMD-check, coverage, and pkgdown; separate queue time from execution time.
- Confirm whether `master` is an active branch or a compatibility trigger.

**Checkpoint:** No workflow trigger is changed until the branch policy and exact required check names are recorded. If `dev` is not PR-gated, revise the plan to keep a full remote check on pushes to `dev`.

### Phase 2: Specify the two CI tiers

- Define the fast feature → `dev` gate: Ubuntu, current supported R release, deterministic package checks/tests, and explicit treatment of `skip_on_cran()` simulation tests.
- Define the full `dev` → `main` gate: existing OS/R matrix, full coverage, pkgdown validation, and stable required-check reporting.
- Specify which post-merge events remain: pkgdown deployment on `main`; only keep full test/coverage push runs if direct-push or post-merge validation requirements justify them.

**Checkpoint:** Confirm that the fast tier covers all required correctness checks and the full tier retains platform-specific, coverage, vignette, and site checks.

### Phase 3: Align local preflight documentation

- Correct `scripts/check-local.sh` comments and relevant contributor guidance so the local check is clearly a preflight and not equivalent to GitHub PR CI.
- Document the fast command, full local command, expected skips, log location, and the remote checks that remain authoritative.
- Do not add an automatic pre-commit hook unless separately requested; keep local execution explicit and predictable.

**Verification:** Run both documented command forms on a clean worktree or verify existing command behavior without changing package code.

### Phase 4: Implement the workflow tiers

- Update workflow conditions/matrices so feature → `dev` gets the fast remote gate and `dev` → `main` gets the full gate.
- Ensure stable required status contexts for both PR targets; use an always-reported aggregate if conditional jobs would otherwise make required checks ambiguous.
- Remove full check/coverage push duplicates only where Phase 1 confirms PR-only branch writes. Keep `main` pkgdown deployment and required release/manual triggers.

**Verification:** Inspect workflow syntax and event/condition logic; confirm paths cannot silently skip an intended required check; run local actionlint if installed and the repository's workflow checks if present.

### Phase 5: Validate, observe, and document

- After explicit user approval, open a small validation PR to `dev`; confirm the fast gate reports its stable required status and completes substantially sooner than the current full matrix. Opening the PR does not authorize merging it.
- Validate the `dev` → `main` route and confirm every full platform/coverage/site gate runs and required checks are enforceable.
- Compare runner minutes and wall time against Phase 1; document the measured change and any tradeoffs in `CLAUDE.md`/CI guidance.

**Completion gate:** All required checks pass on both PR types, no required status is left pending/skipped, branch protection still blocks unverified merges, and the working tree/status reflect the shipped policy. PR creation and merge each require explicit authorization; passing checks alone authorizes neither.

## Risks and mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| `dev` can be updated without a PR | Untested changes enter shared integration | Keep full checks on `dev` pushes until branch protection is confirmed and corrected |
| Conditional jobs change GitHub status names | Merge becomes blocked or a check is no longer enforced | Record exact required contexts first; preserve a stable aggregate status |
| Fast tier skips a defect caught only by another platform or full simulation | Regression reaches `dev` | Retain full matrix/coverage on `dev` → `main`; keep fast tier limited to quick checks |
| Local environment differs from GitHub | Local pass gives false confidence | State that local checks are advisory; remote required checks remain authoritative |
| Docs-only path filters skip a required check | PR remains pending or merges without intended validation | Exclude path filters from first implementation; later use an always-reported aggregate |
| `master` triggers are removed while branch remains active | Branch loses CI coverage | Verify remote branch use before removing compatibility triggers |

## Open questions to resolve in Phase 1

- Which checks are currently required on `dev` and `main`, and can the GitHub integration or maintainer settings view confirm them?
- Are direct pushes, merge queue, or non-PR release pushes part of the actual branch process?
- What are median and high-percentile runner minutes for the matrix, coverage, and pkgdown jobs across recent runs?
- Is `master` still used, or can its triggers be removed?
- Which fast-tier command best matches the package's current test design without running expensive simulations on every feature PR?
