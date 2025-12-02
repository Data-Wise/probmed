# probmed: Probabilistic Effect Sizes for Causal Mediation Analysis

Compute P_med, a scale-free probabilistic effect size for causal
mediation analysis.

Compute P_med, a scale-free probabilistic effect size for causal
mediation analysis. P_med = P(Y\_{x\*,M_x} \> Y\_{x,M_x}) represents the
probability that the outcome under treatment with the mediator at its
control level exceeds the outcome under control.

The `probmed` package provides a robust framework for computing
\\P\_{med}\\, a scale-free probabilistic effect size for causal
mediation analysis. Unlike traditional effect sizes, \\P\_{med}\\
quantifies the probability that the counterfactual outcome under
treatment with the mediator at its control level exceeds the
counterfactual outcome under control.

## Key Features

- **Scale-free interpretation**: \\P\_{med} \in \[0, 1\]\\ has clear
  probabilistic meaning regardless of variable scales

- **Flexible estimation**: Works with GLMs, SEM (via lavaan), and other
  fitted models

- **Comprehensive reporting**: Provides both \\P\_{med}\\ and
  traditional Indirect Effect (\\a \times b\\) with bootstrap confidence
  intervals

- **Multiple mediators**: Supports parallel mediation structures

- **Bootstrap inference**: Parametric and nonparametric bootstrap for
  confidence intervals

- **S7 OOP architecture**: Modern, type-safe object system for
  extensibility

## Main Functions

- [`pmed`](https://data-wise.github.io/probmed/reference/pmed.md):
  Compute \\P\_{med}\\ from formulas or extracted objects

- [`extract_mediation`](https://data-wise.github.io/probmed/reference/extract_mediation.md):
  Extract mediation structure from fitted models

## Supported Model Types

- **lm/glm**: Native support for linear and generalized linear models

- **lavaan**: Structural Equation Models with FIML and robust estimators

- **mediation**: Direct integration with mediate() objects

## Getting Started

The basic workflow involves three steps:

1.  Specify your mediation model (outcome and mediator formulas)

2.  Choose an inference method (parametric bootstrap recommended)

3.  Interpret results (both \\P\_{med}\\ and Indirect Effect)

    # Example workflow
    library(probmed)

    # Generate example data
    set.seed(123)
    n <- 200
    data <- data.frame(X = rnorm(n), C = rnorm(n))
    data$M <- 0.5 * data$X + 0.3 * data$C + rnorm(n)
    data$Y <- 0.4 * data$M + 0.2 * data$X + 0.2 * data$C + rnorm(n)

    # Compute P_med with bootstrap CI
    result <- pmed(
      Y ~ X + M + C,
      formula_m = M ~ X + C,
      data = data,
      treatment = "X",
      mediator = "M",
      method = "parametric_bootstrap",
      n_boot = 1000
    )

    print(result)

## Learn More

- [Package Website](https://data-wise.github.io/probmed/)

- [`vignette("introduction", package = "probmed")`](https://data-wise.github.io/probmed/articles/introduction.md)

- [`vignette("lavaan-integration", package = "probmed")`](https://data-wise.github.io/probmed/articles/lavaan-integration.md)

- [`vignette("mediation-integration", package = "probmed")`](https://data-wise.github.io/probmed/articles/mediation-integration.md)

## References

Tofighi, D., et al. (In Press). Probabilistic Effect Sizes for Mediation
Analysis.

Imai, K., Keele, L., & Tingley, D. (2010). A general approach to causal
mediation analysis. *Psychological Methods*, 15(4), 309-334.

MacKinnon, D. P. (2008). *Introduction to Statistical Mediation
Analysis*. Lawrence Erlbaum Associates.

## See also

Useful links:

- <https://data-wise.github.io/probmed/>

- <https://github.com/data-wise/probmed>

- Report bugs at <https://github.com/data-wise/probmed/issues>

Useful links:

- <https://data-wise.github.io/probmed/>

- <https://github.com/data-wise/probmed>

- Report bugs at <https://github.com/data-wise/probmed/issues>

Useful links:

- <https://data-wise.github.io/probmed/>

- <https://github.com/data-wise/probmed>

- Report bugs at <https://github.com/data-wise/probmed/issues>

## Author

**Maintainer**: Davood Tofighi <dtofighi@gmail.com>
([ORCID](https://orcid.org/0000-0001-8523-7776))
