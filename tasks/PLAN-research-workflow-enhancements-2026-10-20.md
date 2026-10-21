# Plan: Public software claims and numerical evidence identity

Filed 2026-10-20 after explicit user authorization to file the
research-workflow plans.
Status: proposed follow-up; filing only, not implementation or acceptance
evidence.

**Owner:** probmed package maintainer; pmed owns manuscript wording.
**Observed local identity:** `dev` at `a692a43b97a4`. This is a dated
checkout observation, not a release/installed-package assertion; recheck
before execution.

## Scope and holds

Package on hold; preserve I6 stop before I7 and all release/submission holds.
Preserve unrelated work, current scientific rules, active plans and
per-instance external-action gates.
No simulation, installation, manuscript revision, submission, release or
deployment ran for this filing.
This dated plan supplements existing plans; it does not replace their queue or
certify completion.

## Ordered tasks

### 1. Pin the existing parity evidence

- [ ] Reuse the pmed feature-parity report and record checkout SHA,
      DESCRIPTION, installed version/library if inspected, and release/site
      evidence with dates. Read the existing I6 gate before proposing more
      work.
- **Acceptance:** Source, installed, and released identities remain separate;
  no duplicate I6 run or package/manuscript MBCO equivalence is presumed.
- **Verification:** Compare R/support-guard.R, exports/help, pmed report, and
  docs/reports/I6-GATE-2026-10-19.md when available.
- **Dependency:** None; reconcile existing evidence first.
- **Scope:** one evidence record or bounded correction; split implementation
  if it exceeds five files.

### 2. Scope documentation corrections

- [ ] Plan correction of Poisson/count claims and treated-person versus
      control-person wording identified by the existing audit. Disambiguate
      pmed_sensitivity() ratio behavior from manuscript P_med and A2/A5
      sensitivity.
- **Acceptance:** Website/help claims match supported Gaussian/Bernoulli
  additive models; ordinal/count and manuscript-local methods are not
  advertised as implemented.
- **Verification:** Audit current local docs and rendered pkgdown after a
  later authorized change; verify any eventual live site separately.
- **Dependency:** Task 1; retain the existing project approval gates.
- **Scope:** one evidence record or bounded correction; split implementation
  if it exceeds five files.

### 3. Preserve numerical stop and release gates

- [ ] Record which I6 reports/decision files are tracked and which still need
      integration. Keep s6/N=1000 failure under its prespecified rule; do not
      advance I7 from a different run or manuscript gate.
- **Acceptance:** Evidence is traceable to the frozen run; incomplete results
  cannot pass; no hold is lifted by planning.
- **Verification:** Use the owning gate/report, exact IDs and source
  identities; later code changes require the package scientific tests and
  check suite.
- **Dependency:** Task 2; retain the existing project approval gates.
- **Scope:** one evidence record or bounded correction; split implementation
  if it exceeds five files.

## Checkpoint and first action

- [ ] After tasks 1–2, reconcile evidence and unresolved decisions before
      beginning task 3.
- [ ] Record actual verification and remaining owner decisions; do not tick
      future execution from this filing.

**First action:** the named maintainer reads the existing evidence in task 1
and reconciles its
residual scope with the current owning plan. Stop if evidence conflicts or a
scientific/author
decision is unresolved; name that decision rather than silently changing the
plan.

## Coordination

- [Shared Savant enhancement plan][ref1] owns reusable command/skill work.
- [pmed companion plan][ref2] owns its repository changes.
- Future execution must follow this repository's instructions, tests, artifact
  checks and branch model.

[ref1]: </private/tmp/savant-simulation-completeness/tasks/plan-research-workflow-2026-10-20.md>
[ref2]: </Users/davood/projects/research/pmed/docs/PLAN-research-workflow-enhancements-2026-10-20.md>
