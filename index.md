# probmed

PR validation uses the [full local
gate](https://github.com/Data-Wise/probmed/blob/dev/tasks/PLAN-local-ci-transition-2026-10-20.md),
which posts a commit-specific required status after passing. GitHub
Actions still publishes the website on `main`; R-CMD-check remains
required on `main` PRs until branch protection is switched, while
coverage is manual.

**probmed** provides a robust framework for computing **$`P_{med}`$**, a
scale-free probabilistic effect size for causal mediation analysis.

## What is $`P_{med}`$?

$`P_{med}`$ offers an intuitive alternative to traditional mediation
effect sizes. With treatment held at its treated level $`x`$ for both
potential outcomes, it is the probability that the outcome when the
mediator takes its treated level $`M(x)`$ exceeds the outcome when it
takes its control level $`M(x^*)`$, with ties counted as one half:
$`P(Y(x, M(x)) > Y(x, M(x^*))) + \tfrac12 P(Y(x, M(x)) = Y(x, M(x^*)))`$.
The direct effect cancels, so $`P_{med}`$ is mediation-only.

Unlike the “Proportion Mediated” or simple indirect effect coefficients,
$`P_{med}`$: \* **Is Scale-Free**: Interpretable regardless of the units
of measurement for $`X`$, $`M`$, or $`Y`$. \* **Is Bounded**: Always
falls between 0 and 1, with 0.5 meaning no mediation effect, making it
comparable across studies. \* **Handles Non-Linearity**: Can be extended
to Generalized Linear Models (GLMs) where coefficients are not directly
comparable.

## Features

Single-mediator plug-in calls honor an explicit `seed`. MBCO is
deterministic and reports numerical statuses in `result@diagnostics`;
only verified likelihood-ratio endpoints set `converged = TRUE`, and
detected disconnected sets or failed profiles return missing interval
bounds. Parallel mediation keeps its legacy inversion path.

- **Flexible Modeling**: Gaussian (identity link) and binary (logit or
  probit link) mediators and outcomes, with covariates. Binary models
  are standardized over the observed covariate rows. Other families
  (Poisson, ordinal, …), interactions and transformed treatment/mediator
  terms are rejected with an informative error rather than simulated
  incorrectly.
- **Comprehensive Effect Reporting**: Provides both $`P_{med}`$
  (probabilistic effect size) and the traditional Indirect Effect
  ($`a \times b`$) with bootstrap confidence intervals.
- **Robust Inference**: Implements multiple methods for constructing
  confidence intervals:
  - **Parametric Bootstrap**: Efficient and accurate for large samples
    with known distributions.
  - **Nonparametric Bootstrap**: Robust to distributional assumptions
    (resampling-based).
  - **Plug-in Estimator**: Fast point estimation.
  - **MBCO** (Gaussian): deterministic likelihood-ratio interval with
    verified endpoints. It never undercovers but is conservative near
    the double null, flagged by `@diagnostics$pmed$near_double_null`.
  - **Model-based bootstrap** (opt-in, `method = "model_bootstrap"`):
    refits both equations on data simulated from the fitted models, so
    residual scales are re-estimated too. Coverage varies by design; see
    `docs/reports/CALIBRATION-2026-10-14.md` in the repository.
- **Modern Architecture**: Built on the **S7** object-oriented system,
  ensuring type safety, stability, and easy extensibility.
- **Ecosystem Integration**: Works seamlessly with popular R packages:
  - **lavaan**: Extract from Structural Equation Models with FIML and
    robust estimators
  - **mediation**: Direct integration with
    [`mediate()`](https://rdrr.io/pkg/mediation/man/mediate.html)
    objects
  - Standard **lm/glm**: Native support for regression models
- **User-Friendly Interface**: Designed to work with standard R
  `formula` interfaces and `data.frame` inputs.

## Modern estimands

Beyond the core
[`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md), the
package ships cross-fitted, influence-function-based estimators, each
with its own S7 result class. See the [Modern estimands
article](https://data-wise.github.io/probmed/articles/modern-estimands.html)
for details and caveats.

| Function | Result class | What it estimates |
|----|----|----|
| [`ward_residual()`](https://data-wise.github.io/probmed/reference/ward_residual.md) | `GaugePmedResult` | Gauge-calibrated P_med (residual over observed effect) with weak-identification gates |
| [`incr_pmed()`](https://data-wise.github.io/probmed/reference/incr_pmed.md), [`incr_sensitivity()`](https://data-wise.github.io/probmed/reference/incr_sensitivity.md) | `IncrPmedResult` | Incremental P_med and its per-delta tipping analysis |
| [`sobol_pmed()`](https://data-wise.github.io/probmed/reference/sobol_pmed.md) | `SobolPmedResult` | Sobol variance-scale P_med |
| [`wasserstein_pmed()`](https://data-wise.github.io/probmed/reference/wasserstein_pmed.md) | `WassersteinPmedResult` | Wasserstein proportion mediated |
| [`rg_flow_contrast()`](https://data-wise.github.io/probmed/reference/rg_flow_contrast.md) | `RgFlowResult` | Scale-indexed NIE/TE with cluster-robust inference |
| [`pmed_sensitivity()`](https://data-wise.github.io/probmed/reference/pmed_sensitivity.md) | `PmedSensitivity` | Unmeasured mediator-outcome confounding (tipping point) |

## Learn more

- [Cookbook](https://data-wise.github.io/probmed/articles/cookbook.html):
  short recipes for common tasks
- [Choosing a confidence
  interval](https://data-wise.github.io/probmed/articles/confidence-intervals.html):
  what each `method` supports, from the calibration study
- [Sensitivity
  analysis](https://data-wise.github.io/probmed/articles/sensitivity-analysis.html):
  tipping-point analysis for M-Y confounding
- Technical notes: [the P_med
  estimand](https://data-wise.github.io/probmed/articles/technical-note-estimand.html)
  and the [validation
  record](https://data-wise.github.io/probmed/articles/technical-note-validation.html)

## Status and known limitations

Version `0.4.0` (GitHub release only; not submitted to CRAN). For the
binary Wasserstein proportion mediated, estimation near an exactly zero
mediation path is hard: the plug-in bias is of order $`n^{-1/4}`$, none
of the point estimators tested met the frozen bias gate (\|bias\|/SD \<
.10) there, and that gate was shown to be infeasible for any estimator
with values in \[0, 1\] in the zero-direct and weak-direct cells at n =
500, and the Fieller interval is the supported uncertainty statement.
See the [boundary article
section](https://data-wise.github.io/probmed/articles/modern-estimands.html#boundary-behavior)
and `docs/reports/REPORT-boundary-point-risk-2026-10-21.md`.

## Mediationverse Ecosystem

**probmed** is part of the **mediationverse** ecosystem for mediation
analysis in R:

| Package | Purpose | Role |
|----|----|----|
| [**medfit**](https://github.com/data-wise/medfit) | Model fitting, extraction, bootstrap | Foundation |
| **probmed** (this) | Probabilistic effect size (P_med) | Application |
| [**RMediation**](https://github.com/data-wise/rmediation) | Confidence intervals (DOP, MBCO) | Application |
| [**medrobust**](https://github.com/data-wise/medrobust) | Sensitivity analysis | Application |
| [**medsim**](https://github.com/data-wise/medsim) | Simulation infrastructure | Support |

See [Ecosystem
Coordination](https://github.com/data-wise/medfit/blob/main/planning/ECOSYSTEM.md)
for version compatibility and development guidelines.

## Installation

Install from the [Data-Wise
r-universe](https://data-wise.r-universe.dev) — no compiler needed, and
dependencies (e.g. `medfit`) resolve automatically:

``` r

install.packages("probmed", repos = "https://data-wise.r-universe.dev")
```

Or install the development version from
[GitHub](https://github.com/data-wise/probmed):

``` r

# install.packages("devtools")
devtools::install_github("data-wise/probmed")
```

## Quick Example

Here is a basic example showing how to simulate data and compute
$`P_{med}`$ using a parametric bootstrap approach.

### 1. Simulate Data

First, we generate a dataset with a continuous treatment ($`X`$),
mediator ($`M`$), and outcome ($`Y`$), along with a covariate ($`C`$).

``` r

library(probmed)

set.seed(123)
n <- 300
data <- data.frame(X = rnorm(n), C = rnorm(n))

# Mediator model: M ~ X + C
data$M <- 0.5 * data$X + 0.3 * data$C + rnorm(n)

# Outcome model: Y ~ M + X + C
data$Y <- 0.4 * data$M + 0.2 * data$X + 0.2 * data$C + rnorm(n)
```

### 2. Compute $`P_{med}`$

We use the
[`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md)
function to estimate the effect size. We specify the models for the
outcome and mediator, and select the bootstrap method.

``` r

result <- pmed(
  Y ~ X + M + C,             # Outcome model formula
  formula_m = M ~ X + C,     # Mediator model formula
  data = data,
  treatment = "X",
  mediator = "M",
  method = "parametric_bootstrap",
  n_boot = 1000,
  seed = 42
)
```

### 3. View Results

The results object provides a clear summary of both $`P_{med}`$ and the
Indirect Effect with confidence intervals.

``` r

print(result)
#>
#> P_med: Probability of Mediated Shift
#> ====================================
#>
#> Estimate: 0.55
#> 95% CI: [0.53, 0.572]
#>
#> Indirect Effect (Product of Coefficients):
#> Estimate: 0.2
#> 95% CI: [0.129, 0.283]
#>
#> Inference: parametric_bootstrap
#> Bootstrap samples: 1000
#>
#> Treatment contrast: X = 1 vs. X* = 0
#>
#> Interpretation:
#>   P(Y(1, M(1)) > Y(1, M(0))) = 0.55
#>   P that the mediator shift (M(0) -> M(1)) leaves a random individual better off, holding X = 1.

summary(result)
```

**Interpretation**: - **$`P_{med}`$ = 0.55**: Holding treatment at
$`X = 1`$, there is a 55% probability that the outcome under the treated
mediator level exceeds the outcome under the control mediator level.
Since 0.50 represents “no mediation effect”, this indicates a modest
positive mediation effect. - **Indirect Effect = 0.2**: The traditional
product-of-coefficients indirect effect ($`a \times b`$) quantifies the
average change in outcome due to the mediated path.

## Workflow Options

`probmed` offers three workflow approaches to suit your needs:

### 1. Formula Interface (Recommended for Most Users)

Directly specify your models using R formulas:

``` r

result <- pmed(
  Y ~ X + M + C,             # Outcome model
  formula_m = M ~ X + C,     # Mediator model
  data = data,
  treatment = "X",
  mediator = "M",
  method = "parametric_bootstrap"
)
```

### 2. Integration with lavaan (For SEM Users)

Extract from fitted Structural Equation Models with full support for
FIML, robust estimators, and standardized estimates:

``` r

library(lavaan)

# Define SEM model
model <- '
  M ~ a*X + C
  Y ~ b*M + cp*X + C
'

# Fit with FIML for missing data
fit <- sem(model, data = data, missing = "fiml")

# Extract and compute P_med
extract <- extract_mediation(fit, treatment = "X", mediator = "M")
result <- pmed(extract, method = "parametric_bootstrap", n_boot = 5000)
print(result)
```

### 3. Integration with mediation Package

Seamlessly work with objects from the `mediation` package:

``` r

library(mediation)

# Fit models
model_m <- lm(M ~ X + C, data = data)
model_y <- lm(Y ~ M + X + C, data = data)

# Run mediate
med_out <- mediate(model_m, model_y, treat = "X", mediator = "M", boot = TRUE)

# Compute P_med from mediate object
result <- pmed(extract_mediation(med_out))
print(result)
```

## Advanced Features

### Generalized Linear Models (GLMs)

`probmed` supports GLMs for non-normal outcomes. For example, with a
binary outcome:

``` r

# Binary outcome example
data_bin <- data
data_bin$Y_binary <- rbinom(n, 1, plogis(0.4 * data$M + 0.2 * data$X))

result_bin <- pmed(
  Y_binary ~ X + M + C,
  formula_m = M ~ X + C,
  data = data_bin,
  family_y = binomial(),    # Logistic regression for binary outcome
  treatment = "X",
  mediator = "M",
  method = "parametric_bootstrap",
  n_boot = 1000
)
print(result_bin)
```

### Bootstrap Methods Comparison

**When to use each method:**

- **Plugin** (`method = "plugin"`): Fast point estimates only, no
  confidence intervals. Use for quick exploration.
- **Parametric Bootstrap** (`method = "parametric_bootstrap"`): Default
  choice. Fast and accurate when parameter estimates are approximately
  normally distributed.
- **Nonparametric Bootstrap** (`method = "nonparametric_bootstrap"`):
  More robust to distributional assumptions but computationally
  intensive. Use when parametric assumptions are questionable.
- **MBCO** (`method = "mbco"`): Gaussian models only. Deterministic
  likelihood-ratio interval; check `result@diagnostics` for endpoint
  status.
- **Model-based bootstrap** (`method = "model_bootstrap"`, opt-in):
  simulates complete datasets from the fitted models with the observed
  treatment and covariates fixed, refits both equations, and reports the
  percentile interval only when at least 99% of draws are valid. Its
  coverage has not yet been calibrated.

``` r

# Compare methods
result_plugin <- pmed(Y ~ X + M, formula_m = M ~ X, data = data,
                      treatment = "X", mediator = "M", method = "plugin")

result_param <- pmed(Y ~ X + M, formula_m = M ~ X, data = data,
                     treatment = "X", mediator = "M",
                     method = "parametric_bootstrap", n_boot = 1000)

result_nonparam <- pmed(Y ~ X + M, formula_m = M ~ X, data = data,
                        treatment = "X", mediator = "M",
                        method = "nonparametric_bootstrap", n_boot = 1000)
```

## References

- Tofighi, D. (In Press). *Probabilistic Effect Sizes for Mediation
  Analysis*.
- MacKinnon, D. P. (2008). *Introduction to Statistical Mediation
  Analysis*. Lawrence Erlbaum Associates.

## Contributing

We welcome contributions! Please see the [Issue
Tracker](https://github.com/data-wise/probmed/issues) for potential
features or bug reports.
