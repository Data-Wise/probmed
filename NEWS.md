# probmed 0.1.0.9000 (development)

## Features

* `pmed(..., method = "mbco")` adds a deterministic Model-Based Constrained
  Optimization interval (Tofighi & Kelley, 2020): a likelihood-ratio interval
  for both P_med and the indirect effect `a * b`, obtained by inverting the
  constrained-likelihood test rather than by resampling. Gaussian outcome and
  mediator (with covariates) and any contrast `x_ref != x_value`; seed-free and
  grid-resolution-independent. Non-Gaussian models still use the bootstrap
  methods.

## Fixes

* `print(PmedResult)` interpretation line now shows the mediation estimand
  `P(Y(1, M(1)) > Y(1, M(0)))` (manuscript Definition 1) instead of the stale
  direct-effect notation `P(Y_{X*, M_X} > Y_{X, M_X})`. The computed value was
  already correct; only the displayed notation was wrong. A plain-language gloss
  was added beneath it. README cached output updated to match.

# probmed 0.1.0 (2026-06-06)

First GitHub release (non-CRAN).

## Features

* `pmed()` computes P_med — a scale-free probabilistic effect size for causal
  mediation — from a formula or a `medfit::MediationData` object, with plugin,
  parametric-bootstrap, and nonparametric-bootstrap methods, alongside the
  indirect effect (`a * b`).

## Fixes

* `pmed()` now computes the **mediation** estimand
  `P(Y(x, M(x)) > Y(x, M(x*))) + 0.5 P(=)` (manuscript Definition 1): treatment
  held fixed, mediator varied between its treated and control levels with the
  tie term. Previously it computed a direct-effect contrast that depended on the
  direct effect `c'` and could land on the wrong side of 0.5. Verified against
  the closed form `Phi(ab / sqrt(2 (b^2 sigma_M^2 + sigma_Y^2)))` and the
  manuscript `memory_exp` example (`P_med = 0.68`).
* Binary/non-Gaussian outcomes now draw Bernoulli responses through the link
  (previously degenerate, returning 0/1), using the new medfit family slot.
* Parametric bootstrap indexes coefficients **by name** (was selecting
  intercepts as `a`/`b`); nonparametric bootstrap refits on the **correct
  family** (was always Gaussian).

## Ecosystem

* Builds on **medfit (>= 0.3.0)** for model extraction and the family/link slot.
  Part of the mediationverse ecosystem.

---

# probmed 0.0.0.9000

## Major Changes

* **S7 Architecture**: The package has been refactored to use the **S7** object-oriented system for robust class definitions and method dispatch.
* **Website Redesign**: The documentation website now uses the `litera` theme with a clean, academic design matching the `rmediation` package style.
* **Quarto Integration**: The package now uses Quarto (`.qmd`) for the README and vignettes, providing modern publishing capabilities.

## New Features

* **$P_{med}$ Calculation**: Implemented `pmed()` function to compute the probabilistic effect size $P_{med}$ for mediation analysis.
* **GLM Support**: Added support for Generalized Linear Models (e.g., logistic regression) for both mediator and outcome.
* **Bootstrap Inference**: Added parametric and nonparametric bootstrap methods for confidence intervals.
* **lavaan Integration**: Added `extract_mediation()` support for SEM models fitted with the `lavaan` package, including FIML and robust estimators.
* **mediation Integration**: Added support for extracting mediation structures directly from `mediation::mediate()` objects.
* **Indirect Effect Reporting**: `pmed()` now reports the Indirect Effect (product of coefficients) alongside $P_{med}$, including bootstrap confidence intervals.

## Documentation

* **Expanded README**: Added detailed explanation of $P_{med}$, features, and examples.
* **New Vignette**: Added "Introduction to probmed" vignette demonstrating linear and binary outcome examples.
* **Integration Vignettes**: Added dedicated vignettes for `lavaan` and `mediation` package integrations.
* **Comparison Vignette**: Added "Comparing probmed Workflows" to guide users on choosing the best integration method.
