# Implementation Plan: mediationverse compatibility check (2026-10-21)

Existing `tasks/plan.md` / `tasks/todo.md` are historical and still hold unchecked tasks, so this plan uses its own files (repo convention: `plan-<topic>.md`).

## Overview
After the v0.4.0 release, confirm probmed still works against the current medfit and that the ecosystem floors in `DESCRIPTION` and docs are honest. Read-only toward sibling repos (no cross-repo writes without a per-instance request).

## Findings that motivate it
- `DESCRIPTION` requires `medfit (>= 0.3.0)`; installed medfit is 0.7.0; the sibling checkout is 0.10.0 (CRAN version unverified).
- medfit blocks probmed, mediationverse and missingmed in the dependency graph.
- medsim's last commit is a breaking MBCO multiple-imputation default change; probmed references medsim only in `inst/WORDLIST`.
- medfit#213 (`fit_sem` bounds) is open; a probmed memory note relies on it.

## Architecture decisions
- Test probmed against each medfit version in an isolated library, never by modifying the sibling repos.
- A floor change in `DESCRIPTION` is a code change: feature branch, PR, full suite.

## Task list

### Phase 1: Establish the facts (read-only)
- [ ] T1: Record the medfit version on CRAN, in the r-universe, and in the sibling checkout. Acceptance: one table in this plan's Results section. Verify: `available.packages()` and the checkout's `DESCRIPTION`.
- [ ] T2: List every medfit symbol probmed imports or calls (`NAMESPACE`, `R/`). Acceptance: list with the medfit version that introduced each. Verify: grep plus medfit NEWS.

### Checkpoint: facts
- [ ] Floor and API surface are documented; decide whether T3 is needed.

### Phase 2: Test against newer medfit
- [ ] T3: Install medfit 0.10.0 (from the sibling checkout) into a temporary library; run the full suite there. Acceptance: pass/fail/skip counts versus the 827-pass baseline on 0.7.0. Verify: `R_LIBS=<tmp> Rscript -e 'devtools::test()'`.
- [ ] T4: Triage any failure against the baseline (reproduce on installed 0.7.0 first). Acceptance: each failure classified as probmed bug, medfit regression, or environment.

### Checkpoint: compatibility
- [ ] Counts reported; no code changed yet.

### Phase 3: Act on findings (only if needed)
- [ ] T5: If probmed breaks, fix on `feature/*` and open a PR to `dev`; if the floor is wrong, raise or document it. Acceptance: suite green on both medfit versions.
- [x] T6: Record the verified compatibility scope in `.STATUS` and README.

## Risks and mitigations
| Risk | Impact | Mitigation |
|------|--------|------------|
| Installing medfit 0.10.0 changes the user library | Medium | Use a temporary library via `R_LIBS` |
| A medfit API change breaks `extract_mediation()` use | High | T2 lists the surface first |
| Raising the floor excludes CRAN users | Medium | Prefer documenting; raise only if tests need it |

## Open questions
- Is medfit 0.10.0 released on CRAN, or only in the checkout and r-universe?
- Does any probmed simulation script call medsim?

## Results (2026-10-21)

| medfit | Source | Failed | Passed | Warnings |
|--------|--------|-------:|-------:|---------:|
| 0.3.2 | CRAN | 0 | 827 | 0 |
| 0.10.0 | sibling checkout | 2 | 825 | 1 |
| 0.12.0 | r-universe | 2 | 825 | 1 |

API surface (T2): `extract_mediation()`, `MediationData`, `ParallelMediationData`. medsim appears only in `inst/WORDLIST`.

Failure (one test, `test-support-guard.R:98`, missing-rows case): with medfit >= 0.10 the extract keeps only complete rows, so `.check_supported_data()` never sees the NA and `pmed()` returns an estimate (medfit only warns that the models used different rows). On 0.3.2 probmed errors as intended. The offset case is unaffected.

Proposed fix (T5, feature branch, not started): run the missing-value check on the user's `data` in the formula method, before fitting, so it does not depend on what `extract@data` retains.

T5 done on `feature/guard-missing-rows`: `.check_formula_data_complete()` runs on the caller's data before fitting. Full suite after the fix: 832 passed, 0 failed on medfit 0.3.2 and on 0.12.0 (0.10.0 failed only on the guard test, which now passes: guard file run on all three versions).

T6 closed 2026-11-24. [PR #60](https://github.com/Data-Wise/probmed/pull/60)
merged as `ca9b6bf`, and the `v0.4.1` tag contains the pre-fit check. The
declared medfit floor remains `>= 0.3.0`. The README distinguishes that floor
from the versions covered by the saved tests. No tests or simulations were
rerun for T6.

## 2026-11-25 dependency fix

The v0.4.2 candidate pins [medfit v0.20.0](https://github.com/Data-Wise/medfit/releases/tag/v0.20.0)
in `DESCRIPTION` so GitHub installers and CI use the same source. The
declared R, medfit, and S7 minimums stay unchanged. The S7 wrapper help
no longer reports an object length that varies by S7 version.
The package imports `methods::initialize` for S7's S4 registration.
An installed namespace load with only base packages enabled passed
with warnings treated as errors after adding the import.

From an empty R library, pak resolved medfit at commit `4739e8e` and
installed probmed 0.4.2 with medfit 0.20.0 and S7 1.0.0 on R 4.6.2.
The installed package loaded and `pmed(method = "mbco")` on the built-in
`mtcars` data returned 0.316323 with an interval of [0.193987, 0.499470].
The check confirmed that the random-number state stayed the same.

The candidate's [full GitHub matrix](https://github.com/Data-Wise/probmed/actions/workflows/R-CMD-check.yaml?query=branch%3Acodex%2Frelease-0.4.2)
is the hosted check for this dependency pair. Phase 8 remains NO-GO.
I7, package publishing, and pmed submission remain on hold. No new
research simulations ran.
