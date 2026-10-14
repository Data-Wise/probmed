# probmed repair: second implementation batch

Date: 2026-10-09. Parent checkpoint: 3ac3de3. Branch: codex/probmed-repair-2026-10-09. Worktree: /private/tmp/probmed-repair-2026-10-09.

## Scope and changes

Single-mediator Gaussian MBCO now explores both compatible coefficient signs, checks optimizer convergence with one deterministic retry, and reports failed constrained fits as missing likelihoods rather than artificial statistical rejection. Positive variance coordinates preserve the original signed probability constraint. The indirect-effect profile also explores both signs. Exact-null unions from batch 1 remain intact.

The probability profile analytically profiles b at fixed r = Vy/b²: its stationary equation is n*r*b² + ym*b - yy = 0, where yy and ym come from the residualized outcome design. Each sign has one admissible root, leaving two numerical variance coordinates. Cached residual quadratics agree with direct OLS SSE in the fixtures and avoid repeatedly fitting nuisance regressions during domain scans.

Point estimates now use ML residual variances SSE/n, matching the free and constrained Gaussian likelihoods. This intentionally changes earlier residual-df point estimates. The historical test reference was adjusted to the same convention; its assertion tolerance was retained. The independent full-likelihood reference described below supplies a separate check.

A separate bounded-resolution inversion helper scans the full search domain, includes the point estimate and exact-null anchor, verifies root residuals and crossing signs, and records root, truncation, disconnected-set and numerical-failure diagnostics. Failed profiles or detected disconnected sets produce missing reported intervals; truncated bounds are diagnostic search bounds, with converged FALSE. The additive result diagnostics property and print method expose both estimands' statuses. Parallel MBCO retains its legacy inversion behavior and is outside this batch.

## Evidence and review

Before repair, two mocked optimizer-exhaustion tests failed: both constrained profile routines accepted non-converged optimizers. They now return missing likelihoods. Synthetic fixtures cover ordinary roots, full-domain truncation, failed profiles and disconnected accepted components. Signed constraints are reconstructed with positive variances for both treatment-contrast directions.

The independent oracle evaluates the full Gaussian log-likelihood directly from observation residuals and optimizes all coefficients and log variances with nloptr SLSQP and the original unsquared equality. It does not call package likelihood, profiling or constraint helpers. Its deterministic starts cover compatible signs and variance scales; additional variance-start grids give the same tested maxima. Probability profiles agree within 1e-6 per observation for targets .4/.6 and contrasts 1/-2. Independent indirect-effect profiles provide four further assertions. Independently bracketed probability interval endpoints agree within the fixed 1e-4 tolerance for two adjusted-model fixtures. These fixtures support numerical agreement, not a theorem of global optimality.

The initial root tolerance of 1e-6 could fail the fixed likelihood-residual check; tightening it to at most 1e-9 resolved those failures without loosening endpoint or residual gates. The final focused repair-profile suite passes 84 assertions with no failures, warnings or skips. The earlier combined focused run passed 128 assertions before the four additional indirect-effect assertions.

Self-review covered constraint reconstruction, optimizer failure propagation, interval status semantics, compatibility and documentation. No external review of repaired code has run. The bounded scan can miss disconnected components narrower than its recorded spacing; status ok does not prove global connectedness. A future calibration producer must map diagnostic interval status truncated to schema status failed while retaining endpoint status truncated. The existing coverage reader does not perform that mapping.

## Local gates

- `devtools::test(stop_on_failure=TRUE)`: PASS 621; FAIL/WARN/SKIP all 0; 209.1 seconds.
- Fresh temporary installation succeeded at /private/tmp/probmed-repair-library; the global library was unchanged.
- Installed-package smoke fixture: probability .567, interval [.540, .597]; indirect effect .251, interval [.150, .377]; both statuses ok and probability converged TRUE.
- Manuscript baseline harness: explicit seeded estimates .526/.526; both exact-null likelihoods -580.19497988099, matching the independent union maximum; all three probes pass.
- README.qmd rendered successfully to README.md against the candidate library.
- Roxygen-generated targeted help pages retained; unrelated metadata and generated-file changes discarded. nloptr is a Suggests dependency for the independent test oracle, not the runtime profile.
- Eight targeted help files pass tools::checkRd; diff whitespace gate passes.

## Remaining work and limits

Covariance validation, stochastic-control forwarding, observed-row standardization, supported-model and separation guards, and opt-in response bootstrap remain unfinished. Calibration and the costed benchmark have not run. Full CRAN checks, vignette rendering, external code review, PR/CI and integration remain pending. No default bootstrap change, release, push or merge occurred; submission/release holds remain.

Next bounded batch: validate stochastic controls and covariance extraction before adding standardization and response bootstrap under the approved plan.
