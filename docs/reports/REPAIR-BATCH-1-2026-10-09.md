# probmed repair: first implementation batch

Date: 2026-10-09. Base: origin/dev 5c6d837. Branch: codex/probmed-repair-2026-10-09. Worktree: /private/tmp/probmed-repair-2026-10-09.

## Scope

Fix explicit plug-in seed forwarding and both estimands' Gaussian exact-null union. The helper profiles covariates in both reduced equations with ML SSE/n variances. Exact-null branches are selected only for exact zero, without replacing nearby non-null targets. This is the first bounded repair batch, not completion of the expanded plan.

## Evidence

New testthat fixtures failed before changes: three failures, four passing assertions. Repeated explicit seeds produced .536 versus .550; both null likelihoods gave approximately -683.9 versus an independently fitted maximum -582.6 in the covariate-adjusted a=0 fixture. After changes, all seven focused assertions pass, including a fixture favoring b=0.

The temporary installation at /private/tmp/probmed-repair-library passes all three manuscript-repository baseline probes: seeded estimates .526/.526; P_med and indirect-effect exact-null likelihoods both -580.19497988099, equal to the independent union maximum. The global installed package is unchanged. R 4.6.1; medfit 0.4.0; probmed 0.3.0.9000. These checks establish only the tested defects, not profile globality or interval coverage.

Roxygen documentation was regenerated with installed roxygen2 8.1.0; unrelated generated-file and metadata changes were discarded. README.qmd rendered successfully to README.md using the candidate library. NEWS, function help, the introductory vignette source, and status record the bounded change. Vignette rendering and full CRAN checks remain pending; no package release or default change occurred.

## Remaining work

- Validate covariance extraction and independent profile oracle.
- Repair signed-basin profiling, convergence and interval inversion/status behavior.
- Validate and forward stochastic controls across all dispatches.
- Add observed-row standardization, supported-model guards and separation detection.
- Implement the opt-in paired response bootstrap with completion diagnostics.
- Run package/documentation checks and a costed benchmark; obtain resource approval before calibration.

Approved plan/spec and original handoff manifest are copied under planning/repair. They are implementation inputs; the prior external review fingerprints do not certify repaired code.

## Completed local gates

- `devtools::test(filter="repair-baseline", stop_on_failure=TRUE)`: PASS 7; no failures, warnings, or skips.
- `devtools::test(stop_on_failure=TRUE)`: PASS 537; FAIL/WARN/SKIP all 0; 190.1 seconds. This includes the seven new assertions.
- `R CMD INSTALL --library=/private/tmp/probmed-repair-library .`: successful; no global installation changed.
- Manuscript harness `baseline --library /private/tmp/probmed-repair-library`: all three defect probes PASS; exit 0, explicitly not calibration.
- `quarto render README.qmd --to gfm`: successful against temporary candidate library.
- `git diff --check`: successful.

Full CRAN and vignette gates, independent external code review, PR/CI, integration and calibration remain pending. No release, push, or merge was performed.
