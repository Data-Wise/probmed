# CLAUDE.md for probmed

Guidance for Claude Code in this repo. Current state, research lines and holds live in `.STATUS` (single source).

## What this is

R package computing **P_med**, a scale-free probabilistic effect size for causal mediation (P_med in [0, 1]); S7 classes; GLM/lavaan models via **medfit**; Gaussian and binary outcomes.

For treatment X (treated `x`, control `x*`), mediator M, outcome Y:

- P_med = P(Y(x, M(x)) > Y(x, M(x*))) + ½·P(Y(x, M(x)) = Y(x, M(x*)))
- Treatment is held at `x` for **both** potential outcomes; only the mediator varies (draws independent). The direct effect `c'` cancels, so P_med is mediation-only. > 0.5 = positive mediation; 0.5 = null.
- Implemented in `.pmed_core_simple()` (`R/compute-core.R`); matches manuscript Definition 1.

## Commands

```r
remotes::install_deps(dependencies = TRUE)
devtools::load_all(); devtools::document(); devtools::test()
rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "error")
pkgdown::build_site()
```

Run covr with `NOT_CRAN=true`; otherwise `skip_on_cran()` tests are skipped and coverage reads low.

## Coding standards

- R >= 4.1.0, native pipe `|>`, tidyverse style, S7 for classes.
- Always namespace non-base calls: `package::function()`.
- Naming: functions/args/properties snake_case; internals dot-prefixed (`.pmed_core_simple()`); S7 classes CamelCase (`PmedResult`).
- **Scalar branches:** never `ifelse(scalar, vec, vec)` (returns length 1; it made corner-EIF weights one constant per fold from PR #8 to PR #34). Use `if (z == 1) p1 else 1 - p1`.
- Check an estimator by simulation against exact truth with a positive control, never against its own grid.

## Code map

```text
R/
├── aaa-imports.R, probmed-package.R, zzz.R, utils.R   # imports, docs, .onLoad, helpers
├── classes.R, generics.R, data.R                      # PmedResult, pmed() generic, datasets
├── methods-pmed.R, methods-print.R                    # pmed() methods; print/summary/plot
├── compute-core.R                                     # core P_med
├── compute-bootstrap.R, compute-parallel.R            # parametric + nonparametric bootstrap
├── compute-model-bootstrap.R                          # method = "model_bootstrap"
├── compute-mbco.R, compute-mbco-verified.R, compute-parallel-mbco.R  # MBCO intervals
├── support-guard.R, separation.R, standardize.R       # pre-RNG family guard, LP separation, standardization
├── corner.R                                           # shared corner-EIF engine
├── gauge-pmed.R, incremental-pmed.R, incremental-sensitivity.R
└── sobol-pmed.R, wasserstein-pmed.R, rg-flow.R, sensitivity.R
```

Model extraction lives in **medfit** (`extract_mediation()` re-exported; returns `medfit::MediationData`). probmed's bootstrap is its own (`R/compute-bootstrap.R`), not `medfit::bootstrap_mediation()`. Requires medfit >= 0.3.0 (on CRAN).

| Generic | Returns | Notes |
|---------|---------|-------|
| `pmed()` | `PmedResult` | `estimate`, `ci_lower`, `ci_upper`, `boot_estimates`, `ie_estimate` |
| `ward_residual()` | `GaugePmedResult` | gauge-calibrated P_med |
| `incr_pmed()` | `IncrPmedResult` | incremental P_med^delta |
| `sobol_pmed()` | `SobolPmedResult` | Sobol variance-scale |
| `wasserstein_pmed()` | `WassersteinPmedResult` | Wasserstein proportion mediated |
| `rg_flow_contrast()` | `RgFlowResult` | scale-indexed NIE/TE; cluster-robust EIF |
| `pmed_sensitivity()` | `PmedSensitivity` | M-Y sensitivity |

## Computation methods (`pmed(method = )`)

| Method | Description | Use when |
|--------|-------------|----------|
| `"plugin"` | Point estimate only | Quick checks |
| `"parametric_bootstrap"` | Sample from N(theta, Sigma) | Default, fast |
| `"nonparametric_bootstrap"` | Resample data, refit | Normality questionable |
| `"mbco"` | Deterministic LR-inversion interval (Tofighi & Kelley 2020), nloptr inner solver | Gaussian + covariates. Never undercovers; conservative (up to ~1-α²) near the double null; `@diagnostics$pmed$near_double_null` flags its core. Calibrated cutoffs rejected (`docs/reports/MBCO-NULL-T5-FEASIBILITY-2026-10-14.md`) |
| `"model_bootstrap"` | Simulate M, then Y given M*, refit; >= 99% valid draws | Binary M/Y with covariates; coverage: `docs/reports/CALIBRATION-2026-10-14.md` |

## Testing

- Target > 90% overall, 100% for the P_med core.
- `tests/testthat/`: one `test-<topic>.R` per `R/` area (24 files, 5 helpers). `test-pmed-estimand.R` is the estimand gate (closed form, c'-invariance, IE = a*b).

## Ecosystem

probmed is an application package in the mediationverse (medfit foundation, RMediation CIs, medrobust sensitivity, medsim simulation; all under `github.com/data-wise/`). Cascade checklist: `docs/CASCADE-PLAN-2026-06-25.md`. Each package's own `.STATUS` is the source of truth for its state. Manuscript source: `~/projects/research/pmed/` (target *Psychological Methods*).

## Git workflow

- `main` protected (PR required, no force-push/deletion; required check `ubuntu-latest (release)`; admin override available). `dev` = integration. `gh-pages` = pkgdown. `feature/*` = worktrees.
- **CI policy (2026-10-20):** hosted R-CMD-check runs on PRs to `main` only. Commit-bound local CI (PR #53, `local-ci/probmed-full`) is being staged; do not disable the hosted required check until a real PR head has a verified local status (`tasks/PLAN-local-ci-transition-2026-10-20.md`).
- **Holds (author, 2026-10-20/21):** PRs #54/#55 stay unmerged until the frozen boundary diagnostics pass; passing checks do not authorize release, submission, or deployment.

**Version**: 0.3.0.9000 (dev; release on hold pending manuscript submission) · **Last Updated**: 2026-10-21
