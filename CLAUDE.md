# CLAUDE.md for probmed Package

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

------------------------------------------------------------------------

## About This Package

**probmed** is an R package for computing P_med, a scale-free
probabilistic effect size for causal mediation analysis. P_med
quantifies mediation as the probability that the outcome under treatment
with the mediator at its **treated** level exceeds the outcome with the
mediator at its **control** level (treatment held fixed for both, so the
direct effect cancels). See the P_med Definition section below.

> ⚠️ **In-flight (2026-06-19, see `.STATUS`):**
> `feature/incremental-pmed` promotes
> [`incr_pmed()`](https://data-wise.github.io/probmed/reference/incr_pmed.md)/`IncrPmedResult`
> (Paper 2, P_med^δ; 169/0, CI-clean, PR-ready). **Hold the PR until you
> intend to integrate.** Pre-CRAN TODO: g-score EIF term (SE currently
> treats fitted g(C) as fixed → conservative). Companion manuscript:
> `~/projects/research/pmed-modern/02-incremental-pmed`.

### Core Mission

Provide researchers with a scale-free, interpretable effect size measure
(P_med in \[0, 1\]) for mediation analysis that works consistently
across mixed variable types and avoids limitations of traditional
mean-difference approaches.

### Key Features

- Scale-free interpretation: P_med in \[0, 1\] with clear probabilistic
  meaning
- Flexible estimation: Works with GLMs, SEM (lavaan), and fitted models
  via medfit
- Mixed outcome types: continuous (Gaussian) and binary outcomes,
  simulated on the correct link scale via medfit’s family/link slot
- Bootstrap inference: Parametric and nonparametric methods
- S7 OOP architecture for type safety and extensibility

------------------------------------------------------------------------

## Common Development Commands

``` r

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

------------------------------------------------------------------------

## Coding Standards

### R Version and Style

- **Minimum R version**: 4.1.0 (native pipe `|>` support)
- **OOP Framework**: S7 (modern object system)
- **Style**: tidyverse style guide with native pipe
- **Namespacing**: ALWAYS use explicit `package::function()` for
  non-base functions

### Naming Conventions

| Type | Convention | Examples |
|----|----|----|
| Functions | snake_case | [`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md), [`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md) |
| Internal | dot prefix | `.pmed_core()`, [`.pmed_parametric_boot()`](https://data-wise.github.io/probmed/reference/dot-pmed_parametric_boot.md) |
| S7 Classes | CamelCase | `PmedResult` |
| Properties | snake_case | `@a_path`, `@ci_lower` |
| Arguments | snake_case | `formula_y`, `treatment`, `n_boot`, `ci_level` |

### Code Organization

    R/
    ├── aaa-imports.R           # Package imports
    ├── classes.R               # S7 class definitions (PmedResult)
    ├── generics.R              # S7 generic functions (pmed)
    ├── methods-pmed.R          # pmed() methods (formula + MediationData)
    ├── methods-print.R         # print/summary/plot methods
    ├── compute-core.R          # Core P_med computation (.pmed_core_simple)
    ├── compute-bootstrap.R     # Parametric + nonparametric bootstrap
    ├── probmed-package.R       # Package-level docs
    ├── zzz.R                   # .onLoad (S7 dispatch registration)
    └── utils.R                 # Utility functions

Note: model extraction lives in **medfit**
([`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md)
is re-exported); probmed no longer defines its own extraction methods.

------------------------------------------------------------------------

## Code Architecture

### S7 Classes

| Class | Purpose | Key Properties |
|----|----|----|
| `PmedResult` | P_med computation results | `estimate`, `ci_lower`, `ci_upper`, `boot_estimates`, `ie_estimate` |

The extracted mediation structure is
[`medfit::MediationData`](https://data-wise.github.io/medfit/reference/MediationData.html)
(provides `a_path`, `b_path`, `c_prime`, `vcov`, `sigma_m`/`sigma_y`,
and `family_m`/`family_y`) — not a probmed-owned class.

### S7 Generics

| Generic | Purpose | Returns |
|----|----|----|
| [`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md) | Extract from fitted models (re-exported from medfit) | [`medfit::MediationData`](https://data-wise.github.io/medfit/reference/MediationData.html) |
| [`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md) | Compute P_med | `PmedResult` |

### P_med Definition

For treatment X (treated `x`, control `x*`), mediator M, outcome Y: -
P_med = P(Y(x, M(x)) \> Y(x, M(x*))) + ½·P(Y(x, M(x)) = Y(x, M(x*))) -
Treatment is held at `x` for **both** potential outcomes; only the
mediator varies between its treated `M(x)` and control `M(x*)` levels
(drawn independently). The direct effect `c'` cancels, so P_med is
mediation-only. - Values \> 0.5 indicate positive mediation; 0.5 = null
(no mediation effect) - Implemented in
[`.pmed_core_simple()`](https://data-wise.github.io/probmed/reference/dot-pmed_core_simple.md)
(`R/compute-core.R`), matching the manuscript estimand (Definition 1).

### Computation Methods

| Method | Description | Use When |
|----|----|----|
| `"plugin"` | Point estimate only | Quick checks, no CI needed |
| `"parametric_bootstrap"` | Sample from N(theta, Sigma) | Default, fast |
| `"nonparametric_bootstrap"` | Resample data, refit | Normality questionable |
| `"mbco"` | LR-inversion interval (Tofighi & Kelley 2020); deterministic | Gaussian + covariates; closed-form CI for P_med and a\*b |

------------------------------------------------------------------------

## Testing Strategy

### Coverage Targets

- **Target**: \>90% overall, 100% for P_med computation core
- Test P_med accuracy, bootstrap reproducibility, model extraction, edge
  cases

### Test Organization

    tests/testthat/
    ├── test-pmed-glm.R          # Formula/GLM interface + parametric bootstrap CIs
    ├── test-pmed-ie.R           # Indirect effect (a*b) reporting
    └── test-pmed-estimand.R     # Estimand correctness: closed-form agreement,
                                 #   c'-invariance, binary non-degeneracy, IE = a*b

------------------------------------------------------------------------

## Ecosystem Coordination

probmed is an **application package** in the mediationverse ecosystem.

### Central Planning

Ecosystem coordination managed in `/Users/dt/mediation-planning/`: -
`ECOSYSTEM-COORDINATION.md` - Version matrix, release timeline -
`MONTHLY-CHECKLIST.md` - Health checks

### Related Packages

| Package | Repository | Purpose |
|----|----|----|
| medfit | <https://github.com/data-wise/medfit> | Foundation (model fitting, extraction) |
| RMediation | <https://github.com/data-wise/rmediation> | Confidence intervals (DOP, MBCO) |
| medrobust | <https://github.com/data-wise/medrobust> | Sensitivity analysis |
| medsim | <https://github.com/data-wise/medsim> | Simulation infrastructure |

### Integration with medfit

probmed uses medfit for: - Model extraction via
[`medfit::extract_mediation()`](https://data-wise.github.io/medfit/reference/extract_mediation.html)
→
[`medfit::MediationData`](https://data-wise.github.io/medfit/reference/MediationData.html) -
The outcome/mediator family/link (`@family_y`, `@family_m`), used to
simulate non-Gaussian potential outcomes on the correct scale

probmed’s bootstrap (parametric + nonparametric) is implemented in
`R/compute-bootstrap.R` — it does **not** use
[`medfit::bootstrap_mediation()`](https://data-wise.github.io/medfit/reference/bootstrap_mediation.html).

Requires **medfit (\>= 0.3.0)**; the GitHub build pins
`Remotes: data-wise/medfit@v0.3.0`.

------------------------------------------------------------------------

## Key References

- McGraw & Wong (1992): Common language effect size statistics
- Vargha & Delaney (2000): CL effect sizes improvements
- Imai, Keele, & Tingley (2010): Causal mediation analysis
- Tofighi et al. (2025): P_med manuscript (in preparation)
  - Manuscript source: `~/projects/research/pmed/` (pmed.qmd, apaquarto;
    target *Psychological Methods*)

------------------------------------------------------------------------

## Git Workflow

Multi-branch pattern (aligned with global craft-style workflow):

- `main` - production (protected: PR required, no force-push, no
  deletions; required status check: `ubuntu-latest (release)`)
- `dev` - integration / planning hub
- `gh-pages` - pkgdown site
- `feature/*` - isolated work in worktrees

Branch protection on `main` was added 2026-05-09 via
`gh api repos/Data-Wise/probmed/branches/main/protection`. Admin
override available (`enforce_admins: false`).

------------------------------------------------------------------------

## Development Roadmap

### Completed

S7 class architecture (`PmedResult`)

Formula interface with GLM

Bootstrap inference (parametric + nonparametric)

medfit integration (extraction, family/link slot)

Correct mediation estimand (matches manuscript Definition 1)

Binary / non-Gaussian outcomes (Bernoulli simulation via the link)

Print/summary/plot methods

method = “mbco” — closed-form MBCO interval (P_med + indirect effect)

### Future

P_med for parallel / sequential structures (medfit already provides
`ParallelMediationData` / `SerialMediationData`; probmed needs the
estimands)

Exposure-mediator interactions (`InteractionMediationData`)

Expand vignettes; integration with medrobust

------------------------------------------------------------------------

**Version**: 0.3.0 (release-prep; non-CRAN GitHub release, pins
<medfit@v0.3.0>) · **Last Updated**: 2026-06-21

> Active research lines beyond the merged core (see `.STATUS` for
> worktree detail): gauge-calibrated P_med
> (`ward_residual`/`GaugePmedResult`, PR \#8 + bootstrap-SE follow-up),
> incremental P_med^δ (`incr_pmed`, g-score EIF merged via PR \#10),
> Sobol variance-scale (`sobol_pmed`), Wasserstein, and the shared M–Y
> sensitivity helper (`pmed_sensitivity`).
