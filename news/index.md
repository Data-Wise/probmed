# Changelog

## probmed (development version)

### Ecosystem Notes

- Part of the mediationverse ecosystem for mediation analysis
- Builds on medfit for model infrastructure (integration planned)
- See [Ecosystem
  Coordination](https://github.com/data-wise/medfit/blob/main/planning/ECOSYSTEM.md)
  for guidelines

------------------------------------------------------------------------

## probmed 0.0.0.9000

### Major Changes

- **S7 Architecture**: The package has been refactored to use the **S7**
  object-oriented system for robust class definitions and method
  dispatch.
- **Website Redesign**: The documentation website now uses the `litera`
  theme with a clean, academic design matching the `rmediation` package
  style.
- **Quarto Integration**: The package now uses Quarto (`.qmd`) for the
  README and vignettes, providing modern publishing capabilities.

### New Features

- **$P_{med}$ Calculation**: Implemented
  [`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md)
  function to compute the probabilistic effect size $P_{med}$ for
  mediation analysis.
- **GLM Support**: Added support for Generalized Linear Models (e.g.,
  logistic regression) for both mediator and outcome.
- **Bootstrap Inference**: Added parametric and nonparametric bootstrap
  methods for confidence intervals.
- **lavaan Integration**: Added
  [`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md)
  support for SEM models fitted with the `lavaan` package, including
  FIML and robust estimators.
- **mediation Integration**: Added support for extracting mediation
  structures directly from
  [`mediation::mediate()`](https://rdrr.io/pkg/mediation/man/mediate.html)
  objects.
- **Indirect Effect Reporting**:
  [`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md) now
  reports the Indirect Effect (product of coefficients) alongside
  $P_{med}$, including bootstrap confidence intervals.

### Documentation

- **Expanded README**: Added detailed explanation of $P_{med}$,
  features, and examples.
- **New Vignette**: Added “Introduction to probmed” vignette
  demonstrating linear and binary outcome examples.
- **Integration Vignettes**: Added dedicated vignettes for `lavaan` and
  `mediation` package integrations.
- **Comparison Vignette**: Added “Comparing probmed Workflows” to guide
  users on choosing the best integration method.
