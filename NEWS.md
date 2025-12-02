# probmed 0.0.0.9000

## Major Changes

*   **S7 Architecture**: The package has been refactored to use the **S7** object-oriented system for robust class definitions and method dispatch.
*   **Website Redesign**: The documentation website now uses the `litera` theme with a clean, academic design matching the `rmediation` package style.
*   **Quarto Integration**: The package now uses Quarto (`.qmd`) for the README and vignettes, providing modern publishing capabilities.

## New Features

*   **$P_{med}$ Calculation**: Implemented `pmed()` function to compute the probabilistic effect size $P_{med}$ for mediation analysis.
*   **GLM Support**: Added support for Generalized Linear Models (e.g., logistic regression) for both mediator and outcome.
*   **Bootstrap Inference**: Added parametric and nonparametric bootstrap methods for confidence intervals.

## Documentation

*   **Expanded README**: Added detailed explanation of $P_{med}$, features, and examples.
*   **New Vignette**: Added "Introduction to probmed" vignette demonstrating linear and binary outcome examples.
