# CLAUDE.md for probmed Package

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## About This Package

**probmed** is an R package for computing P_med, a scale-free probabilistic effect size for causal mediation analysis. P_med quantifies mediation as the probability that the outcome under treatment with the mediator at its **treated** level exceeds the outcome with the mediator at its **control** level (treatment held fixed for both, so the direct effect cancels). See the P_med Definition section below.

### Core Mission

Provide researchers with a scale-free, interpretable effect size measure (P_med in [0, 1]) for mediation analysis that works consistently across mixed variable types and avoids limitations of traditional mean-difference approaches.

### Key Features

- Scale-free interpretation: P_med in [0, 1] with clear probabilistic meaning
- Flexible estimation: Works with GLMs, SEM (lavaan), and fitted models via medfit
- Mixed outcome types: continuous (Gaussian) and binary outcomes, simulated on the
  correct link scale via medfit's family/link slot
- Bootstrap inference: Parametric and nonparametric methods
- S7 OOP architecture for type safety and extensibility

---

## Common Development Commands

```r
# Install dependencies and check package
remotes::install_deps(dependencies = TRUE)
rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "error")

# Development workflow
devtools::load_all()
devtools::document()
devtools::test()

# Build pkgdown site
pkgdown::build_site()
```

---

## Coding Standards

### R Version and Style

- **Minimum R version**: 4.1.0 (native pipe `|>` support)
- **OOP Framework**: S7 (modern object system)
- **Style**: tidyverse style guide with native pipe
- **Namespacing**: ALWAYS use explicit `package::function()` for non-base functions
- **Scalar branches**: never `ifelse(scalar, vec, vec)` — it returns length 1 (the
  corner-EIF weights were one constant per fold from PR #8 to PR #34 because of this).
  Use `if (z == 1) p1 else 1 - p1`. Simulation checks of an estimator run against
  exact truth with a positive control, not against the estimator's own grid.

### Naming Conventions

| Type | Convention | Examples |
|------|------------|----------|
| Functions | snake_case | `pmed()`, `extract_mediation()` |
| Internal | dot prefix | `.pmed_core()`, `.pmed_parametric_boot()` |
| S7 Classes | CamelCase | `PmedResult` |
| Properties | snake_case | `@a_path`, `@ci_lower` |
| Arguments | snake_case | `formula_y`, `treatment`, `n_boot`, `ci_level` |

### Code Organization

```
R/
├── aaa-imports.R, probmed-package.R, zzz.R, utils.R   # imports, package docs, .onLoad, helpers
├── classes.R, generics.R, data.R                      # PmedResult, pmed() generic, datasets
├── methods-pmed.R, methods-print.R                    # pmed() methods; print/summary/plot
├── compute-core.R                                     # Core P_med (.pmed_core_simple)
├── compute-bootstrap.R, compute-parallel.R            # Parametric + nonparametric bootstrap
├── compute-model-bootstrap.R                          # method = "model_bootstrap"
├── compute-mbco.R, compute-mbco-verified.R, compute-parallel-mbco.R  # MBCO intervals
├── support-guard.R, separation.R, standardize.R       # Pre-RNG family guard, LP separation, standardization
├── corner.R                                           # Shared corner-EIF engine
├── gauge-pmed.R                                       # ward_residual / GaugePmedResult
├── incremental-pmed.R, incremental-sensitivity.R      # incr_pmed, incr_sensitivity
├── sobol-pmed.R, wasserstein-pmed.R                   # sobol_pmed, wasserstein_pmed
├── rg-flow.R                                          # rg_flow_contrast / RgFlowResult
└── sensitivity.R                                      # pmed_sensitivity / PmedSensitivity
```

Note: model extraction lives in **medfit** (`extract_mediation()` is re-exported);
probmed no longer defines its own extraction methods.

---

## Code Architecture

### S7 Classes

| Class | Purpose | Key Properties |
|-------|---------|----------------|
| `PmedResult` | P_med computation results | `estimate`, `ci_lower`, `ci_upper`, `boot_estimates`, `ie_estimate` |
| `GaugePmedResult` | Gauge-calibrated P_med | `ward_residual()` output |
| `IncrPmedResult` | Incremental P_med^delta | `incr_pmed()` output |
| `SobolPmedResult` | Sobol variance-scale P_med | `sobol_pmed()` output |
| `WassersteinPmedResult` | Wasserstein proportion mediated | `wasserstein_pmed()` output |
| `RgFlowResult` | Scale-indexed proportion mediated (NIE/TE) | `rg_flow_contrast()` output |
| `PmedSensitivity` | M-Y sensitivity | `pmed_sensitivity()` output |

The extracted mediation structure is `medfit::MediationData` (provides `a_path`,
`b_path`, `c_prime`, `vcov`, `sigma_m`/`sigma_y`, and `family_m`/`family_y`) — not
a probmed-owned class.

### S7 Generics

| Generic | Purpose | Returns |
|---------|---------|---------|
| `extract_mediation()` | Extract from fitted models (re-exported from medfit) | `medfit::MediationData` |
| `pmed()` | Compute P_med | `PmedResult` |
| `ward_residual()` | Gauge-calibrated P_med | `GaugePmedResult` |
| `incr_pmed()` | Incremental P_med^delta | `IncrPmedResult` |
| `sobol_pmed()` | Sobol variance-scale P_med | `SobolPmedResult` |
| `wasserstein_pmed()` | Wasserstein proportion mediated | `WassersteinPmedResult` |

### P_med Definition

For treatment X (treated `x`, control `x*`), mediator M, outcome Y:
- P_med = P(Y(x, M(x)) > Y(x, M(x*))) + ½·P(Y(x, M(x)) = Y(x, M(x*)))
- Treatment is held at `x` for **both** potential outcomes; only the mediator
  varies between its treated `M(x)` and control `M(x*)` levels (drawn
  independently). The direct effect `c'` cancels, so P_med is mediation-only.
- Values > 0.5 indicate positive mediation; 0.5 = null (no mediation effect)
- Implemented in `.pmed_core_simple()` (`R/compute-core.R`), matching the
  manuscript estimand (Definition 1).

### Computation Methods

| Method | Description | Use When |
|--------|-------------|----------|
| `"plugin"` | Point estimate only | Quick checks, no CI needed |
| `"parametric_bootstrap"` | Sample from N(theta, Sigma) | Default, fast |
| `"nonparametric_bootstrap"` | Resample data, refit | Normality questionable |
| `"mbco"` | LR-inversion interval (Tofighi & Kelley 2020); deterministic; nloptr inner solver | Gaussian + covariates; closed-form CI for P_med and a*b. Never undercovers; conservative (up to ~1-α²) while both paths are within a few SE of zero — `@diagnostics$pmed$near_double_null` flags the core of that band. Calibrated cutoffs rejected (Phase 8, `docs/reports/MBCO-NULL-T5-FEASIBILITY-2026-10-14.md`) |
| `"model_bootstrap"` | Response bootstrap: simulate M, then Y given M*, refit; ≥99% valid draws | Binary M/Y with covariates; see `docs/reports/CALIBRATION-2026-10-14.md` for where coverage holds |

---

## Testing Strategy

### Coverage Targets

- **Target**: >90% overall, 100% for P_med computation core
- Test P_med accuracy, bootstrap reproducibility, model extraction, edge cases

### Test Organization

```
tests/testthat/   # 24 test files, 5 helper-*.R; one test-<topic>.R per R/ area
├── test-pmed-estimand.R   # Estimand correctness: closed form, c'-invariance, IE = a*b (core, 100%)
├── test-pmed-glm.R, test-pmed-ie.R, test-pmed-print.R, test-pmed-parallel.R
├── test-pmed-mbco.R, test-mbco-null*.R, test-model-bootstrap.R, test-crn.R
├── test-support-guard.R, test-separation.R, test-standardize.R
└── test-gauge-pmed.R, test-incremental-*.R, test-sobol-pmed.R, test-wasserstein-pmed.R,
    test-rg-flow.R, test-sensitivity.R, test-repair-*.R, test-covariance-extraction.R,
    test-data-multilevel-designs.R
```
Run covr with `NOT_CRAN=true`, otherwise `skip_on_cran()` tests are skipped and coverage reads low.

---

## Ecosystem Coordination

probmed is an **application package** in the mediationverse ecosystem.

### Central Planning

Cascade checklist: `docs/CASCADE-PLAN-2026-06-25.md`. Each package's own `.STATUS` is the source of truth for its state.

### Related Packages

| Package | Repository | Purpose |
|---------|-----------|---------|
| medfit | https://github.com/data-wise/medfit | Foundation (model fitting, extraction) |
| RMediation | https://github.com/data-wise/rmediation | Confidence intervals (DOP, MBCO) |
| medrobust | https://github.com/data-wise/medrobust | Sensitivity analysis |
| medsim | https://github.com/data-wise/medsim | Simulation infrastructure |

### Integration with medfit

probmed uses medfit for:
- Model extraction via `medfit::extract_mediation()` → `medfit::MediationData`
- The outcome/mediator family/link (`@family_y`, `@family_m`), used to simulate
  non-Gaussian potential outcomes on the correct scale

probmed's bootstrap (parametric + nonparametric) is implemented in
`R/compute-bootstrap.R` — it does **not** use `medfit::bootstrap_mediation()`.

Requires **medfit (>= 0.3.0)**; medfit 0.3.2 is on CRAN, so the former
`Remotes: data-wise/medfit@v0.3.0` pin was dropped (`b552b88`, 2026-07-24).

---

## Key References

- McGraw & Wong (1992): Common language effect size statistics
- Vargha & Delaney (2000): CL effect sizes improvements
- Imai, Keele, & Tingley (2010): Causal mediation analysis
- Tofighi et al. (2025): P_med manuscript (in preparation)
  - Manuscript source: `~/projects/research/pmed/` (pmed.qmd, apaquarto; target *Psychological Methods*)

---

## Git Workflow

Multi-branch pattern (aligned with global craft-style workflow):

- `main` - production (protected: PR required, no force-push, no deletions; required status check: `ubuntu-latest (release)`)
- `dev` - integration / planning hub
- `gh-pages` - pkgdown site
- `feature/*` - isolated work in worktrees

Branch protection on `main` was added 2026-05-09 via `gh api repos/Data-Wise/probmed/branches/main/protection`. Admin override available (`enforce_admins: false`).

**CI policy (2026-10-20):** hosted R-CMD-check runs on PRs to `main` only. A commit-bound local CI (PR #53, `local-ci/probmed-full`) is being staged; do not disable the hosted required check until a real PR head has a verified local status. Details: `tasks/PLAN-local-ci-transition-2026-10-20.md`.

**Holds (author, 2026-10-20/21):** PRs #54/#55 stay unmerged until the frozen boundary diagnostics pass; no release, submission, or deployment is authorized by passing checks. See `.STATUS`.

---

## Development Roadmap

### Completed

- [x] S7 class architecture (`PmedResult`)
- [x] Formula interface with GLM
- [x] Bootstrap inference (parametric + nonparametric)
- [x] medfit integration (extraction, family/link slot)
- [x] Correct mediation estimand (matches manuscript Definition 1)
- [x] Binary / non-Gaussian outcomes (Bernoulli simulation via the link)
- [x] Print/summary/plot methods
- [x] method = "mbco" — closed-form MBCO interval (P_med + indirect effect)

### Future
- [ ] P_med for parallel / sequential structures (medfit already provides
      `ParallelMediationData` / `SerialMediationData`; probmed needs the estimands)
- [ ] Exposure-mediator interactions (`InteractionMediationData`)
- [ ] Expand vignettes; integration with medrobust

---

**Version**: 0.3.0.9000 (dev; release ON HOLD pending manuscript submission) · **Last Updated**: 2026-10-21

Current state, active research lines and holds: `.STATUS` (single source; not restated here).
