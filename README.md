
<!-- README.md is generated from README.qmd. Please edit that file -->

# probmed

<!-- badges: start -->

[![R-CMD-check](https://github.com/data-wise/probmed/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/data-wise/probmed/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/data-wise/probmed/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/data-wise/probmed/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

The goal of `probmed` is to provide a robust and flexible framework for
computing **P_med**, a scale-free probabilistic effect size for causal
mediation analysis. Unlike traditional effect sizes, P_med offers an
intuitive interpretation of the mediation effect in terms of
probability, making it applicable across various model types (linear,
GLM, etc.).

## Features

- **Scale-Free Effect Size**: Computes P_med, which is interpretable
  regardless of the scale of the variables.
- **Flexible Modeling**: Supports both linear and generalized linear
  models (GLMs) for mediator and outcome.
- **Modern Architecture**: Built using the **S7** object-oriented system
  for robustness and extensibility.
- **Inference Methods**: Includes parametric bootstrap, nonparametric
  bootstrap, and plug-in estimators.

## Installation

You can install the development version of probmed from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dtofighi/probmed")
```

## Quick Example

Here is a basic example showing how to simulate data and compute P_med
using a parametric bootstrap approach:

``` r
library(probmed)

# 1. Simulate data
set.seed(123)
n <- 300
data <- data.frame(X = rnorm(n), C = rnorm(n))
# Mediator model: M ~ X + C
data$M <- 0.5 * data$X + 0.3 * data$C + rnorm(n)
# Outcome model: Y ~ M + X + C
data$Y <- 0.4 * data$M + 0.2 * data$X + 0.2 * data$C + rnorm(n)

# 2. Compute P_med
# We specify the outcome model formula (first arg) and the mediator model formula.
result <- pmed(
  Y ~ X + M + C,
  formula_m = M ~ X + C,
  data = data,
  treatment = "X",
  mediator = "M",
  method = "parametric_bootstrap",
  n_boot = 1000,
  seed = 42
)

# 3. View Results
print(result)
summary(result)
```

## Advanced Usage

`probmed` allows you to specify different families for GLMs (e.g.,
binomial for binary outcomes) and customize the reference values for the
treatment effect.

``` r
# Example with binary outcome (not run)
# result_bin <- pmed(
#   Y ~ X + M + C,
#   formula_m = M ~ X + C,
#   data = data_bin,
#   family_y = binomial(),
#   treatment = "X",
#   mediator = "M"
# )
```
