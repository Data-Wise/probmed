# Plan: Public software claims and numerical evidence identity

Filed 2026-10-20 after explicit user authorization to file the
research-workflow plans.
Status: parity audit complete and package documentation corrections merged
to dev. Verification of the published corrections is pending an authorized
site deployment.

**Owner:** probmed package maintainer; pmed owns manuscript wording.
**Observed local identity:** `dev` at `a692a43b97a4`. This is a dated
checkout observation, not a release/installed-package assertion; recheck
before execution.

## Scope and holds

Phase 8 remains NO-GO. Preserve the I6 stop before I7 and the holds on the
package release and pmed submission.
Preserve unrelated work, current scientific rules, active plans and
per-instance external-action gates.
No simulation, installation, manuscript revision, submission, release or
deployment ran for this filing.
This dated plan supplements existing plans; it does not replace their queue or
certify completion.

## Ordered tasks

### 1. Pin the existing parity evidence

- [x] Reuse the pmed feature-parity report and record checkout SHA,
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

**2026-11-23 audit:** Compared probmed dev `767f623` (DESCRIPTION version
`0.4.1`) with pmed draft `29af925`. The installed package reported `0.4.1`
from `505e331`. R-universe served `0.4.1` from `ab411df`, and the live
pkgdown site displayed `0.4.1`. The manuscript's Stage 1 and the live
homepage already scoped count models correctly. The two vignette
descriptions and the sensitivity helper's scope still needed correction.

### 2. Scope documentation corrections

- [x] Correct the remaining treatment-contrast wording and disambiguate the
      `pmed_sensitivity()` ratio from manuscript P_med and A2/A5 sensitivity.
- **Acceptance:** Website/help claims match supported Gaussian/Bernoulli
  additive models; ordinal/count and manuscript-local methods are not
  advertised as implemented.
- **Verification:** Audit current local docs and rendered pkgdown after a
  later authorized change; verify any eventual live site separately.
- **Dependency:** Task 1; retain the existing project approval gates.
- **Scope:** one evidence record or bounded correction; split implementation
  if it exceeds five files.

**2026-11-23 merge:** [PR #71](https://github.com/Data-Wise/probmed/pull/71)
squash-merged reviewed head `bb37e2a` into dev as `bc16bfd`. The squash tree
matched the reviewed head. The four edited vignettes and both pkgdown
reference pages were rendered and inspected locally. R CMD check reported
zero errors, warnings, and notes with tests, examples, and vignette execution
disabled. No new simulations ran.

### 3. Preserve numerical stop and release gates

- [x] Confirm the [I6 report](../docs/reports/I6-GATE-2026-10-19.md) and
      [decision CSV](../docs/reports/I6-GATE-DECISIONS-2026-10-19.csv) are
      tracked. The prespecified rule still fails at s6/N=1,000. I7 remains
      stopped.
- **Acceptance:** Evidence is traceable to the frozen run; incomplete results
  cannot pass; no hold is lifted by planning.
- **Verification:** Use the owning gate/report, exact IDs and source
  identities; later code changes require the package scientific tests and
  check suite.
- **Dependency:** Task 2; retain the existing project approval gates.
- **Scope:** one evidence record or bounded correction; split implementation
  if it exceeds five files.

## Checkpoint

- [x] Reconcile the claims audit, package documentation fixes, and existing
      I6 decision in [`.STATUS`](../.STATUS).
- [ ] Verify the corrected vignettes and sensitivity help on the public
      pkgdown site after a separately authorized deployment.

## Coordination

- [Shared Savant enhancement plan][ref1] owns reusable command/skill work.
- [pmed companion plan][ref2] owns its repository changes.
- Future execution must follow this repository's instructions, tests, artifact
  checks and branch model.

[ref1]: </private/tmp/savant-simulation-completeness/tasks/plan-research-workflow-2026-10-20.md>
[ref2]: </Users/davood/projects/research/pmed/docs/PLAN-research-workflow-enhancements-2026-10-20.md>
