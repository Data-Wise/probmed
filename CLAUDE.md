# CLAUDE.md for probmed Package

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## About This Package

probmed is an R package for computing P_med, a scale-free probabilistic
effect size for causal mediation analysis. P_med quantifies mediation
effects as the probability that the counterfactual outcome under
treatment with the mediator at its control level exceeds the
counterfactual outcome under control: P(Y\_{x\*,M_x} \> Y\_{x,M_x}).

### Core Mission

Provide researchers with a scale-free, interpretable effect size measure
for mediation analysis that works consistently across mixed variable
types (binary, continuous) and avoids the limitations of traditional
mean-difference approaches. Methods accommodate both simple and complex
mediation models while maintaining accessibility for applied
practitioners.

### Key Features

- **Scale-free interpretation**: P_med ∈ \[0, 1\] has clear
  probabilistic meaning regardless of variable scales
- **Flexible estimation**: Works with GLMs, SEM (via lavaan), and other
  fitted models
- **Multiple mediators**: Supports parallel and sequential mediation
  structures
- **Interactions**: Handles exposure-mediator interactions with
  conditioning
- **Bootstrap inference**: Parametric and nonparametric bootstrap for
  confidence intervals
- **S7 OOP architecture**: Modern, type-safe object system for
  extensibility

### Key References

- McGraw & Wong (1992): Common language effect size statistics
  (*Psychological Bulletin*)
- Vargha & Delaney (2000): Critique and improvement of CL effect sizes
  (*Journal of Educational and Behavioral Statistics*)
- Imai, Keele, & Tingley (2010): General approach to causal mediation
  analysis (*Psychological Methods*)
- VanderWeele & Vansteelandt (2014): Mediation analysis with multiple
  mediators (*European Journal of Epidemiology*)
- Tofighi et al. (2025): P_med manuscript (target: *Psychological
  Methods*)

## Common Development Commands

### Package Building and Checking

``` r
# Install package dependencies
install.packages(c("remotes", "rcmdcheck"))
remotes::install_deps(dependencies = TRUE)

# Check package (standard R CMD check)
rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "error")

# Build package
devtools::build()

# Install and reload package during development
devtools::load_all()
```

### Documentation

``` r
# Generate documentation from roxygen2 comments
devtools::document()

# Build PDF manual
devtools::build_manual()

# Build pkgdown site
pkgdown::build_site()
```

### Testing and Linting

The package uses GitHub Actions for CI/CD with workflows defined in
`.github/workflows/`:

- `R-CMD-check.yaml`: R package check on multiple OS and R versions
- `test-coverage.yaml`: Test coverage via covr
- `pkgdown.yaml`: Documentation website deployment

## Coding Standards

### R Version and Style

- **Minimum R version**: 4.1.0 (native pipe \|\> support)
- **OOP Framework**: S7 (modern object system)
- **Style**: tidyverse style guide with native pipe
- **Roxygen2**: Use roxygen2 for documentation (\>= 7.2.0)
- **Clarity priority**: Code must be readable for both methodologists
  and applied users

### Naming Conventions

The package uses **snake_case** consistently following tidyverse
conventions:

**Functions:**

- Main exports:
  [`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md),
  [`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md)
- S7 generics follow convention:
  [`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md),
  [`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md)
- Internal functions prefix with dot: `.pmed_core()`,
  [`.pmed_parametric_boot()`](https://data-wise.github.io/probmed/reference/dot-pmed_parametric_boot.md)

**Arguments:**

- `formula_y`, `formula_m` for model specifications
- `treatment`, `mediator`, `outcome` for variable names
- `x_ref`, `x_value` for treatment contrast
- `method` for inference method: “parametric_bootstrap”,
  “nonparametric_bootstrap”, “plugin”
- `n_boot` for number of bootstrap samples
- `ci_level` for confidence level (default: 0.95)

**S7 Classes:**

- CamelCase: `MediationExtract`, `PmedResult`
- Properties use snake_case: `@a_path`, `@b_path`, `@ci_lower`

### Code Organization

    R/
    ├── aaa-imports.R           # Package imports and setup
    ├── classes.R               # S7 class definitions
    ├── generics.R              # S7 generic functions
    ├── methods-extract.R       # extract_mediation() methods
    ├── methods-pmed.R          # pmed() methods
    ├── methods-print.R         # print(), summary(), plot() methods
    ├── compute-core.R          # Core P_med computation
    ├── compute-bootstrap.R     # Bootstrap implementations
    └── utils.R                 # Utility functions

## Code Architecture

### S7 Object System

The package uses S7 for type-safe, modern object-oriented programming:

**Key S7 Classes:**

1.  **`MediationExtract`** - Extracted mediation structure from fitted
    models
    - Properties: `estimates`, `a_path`, `b_path`, `c_prime`, `vcov`,
      `data`, etc.
    - Validation rules ensure internal consistency
    - Subclasses: `MediationExtractLavaan` for lavaan-specific
      properties
2.  **`PmedResult`** - P_med computation results
    - Properties: `estimate`, `ci_lower`, `ci_upper`, `method`,
      `boot_estimates`, etc.
    - Validator checks P_med ∈ \[0, 1\] and CI consistency
    - Source tracking via `source_extract` property

**S7 Generics:**

1.  **[`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md)** -
    Extract mediation structure from models
    - Methods: [`stats::lm`](https://rdrr.io/r/stats/lm.html), `lavaan`,
      `mediation`, etc.
    - Returns: `MediationExtract` object
2.  **[`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md)** -
    Compute P_med from various inputs
    - Methods: `formula`, `MediationExtract`
    - Returns: `PmedResult` object

### Core Function Hierarchy

**User-Facing Functions:**

1.  **`pmed(formula_y, formula_m, ...)`** - Formula interface (most
    common)
    - Fits GLMs automatically
    - Dispatches to appropriate computation method
    - Returns `PmedResult` with CI
2.  **`pmed(object, ...)`** - Works with `MediationExtract`
    - For pre-extracted structures
    - More flexible for custom models
3.  **`extract_mediation(object, ...)`** - Extract from fitted models
    - Generic with methods for lm/glm, lavaan, mediation, medflex
    - Standardizes across packages
    - Returns `MediationExtract`

**Internal Computation Functions:**

- [`.pmed_compute()`](https://data-wise.github.io/probmed/reference/dot-pmed_compute.md):
  Main dispatcher
- [`.pmed_plugin()`](https://data-wise.github.io/probmed/reference/dot-pmed_plugin.md):
  Point estimate only
- [`.pmed_parametric_boot()`](https://data-wise.github.io/probmed/reference/dot-pmed_parametric_boot.md):
  Parametric bootstrap
- [`.pmed_nonparametric_boot()`](https://data-wise.github.io/probmed/reference/dot-pmed_nonparametric_boot.md):
  Nonparametric bootstrap
- [`.pmed_core_simple()`](https://data-wise.github.io/probmed/reference/dot-pmed_core_simple.md):
  Core P_med calculation

### Key Dependencies

**Required:**

- **S7**: Modern object system
- **stats**: GLM fitting and prediction
- **methods**: S3/S4 compatibility checks
- **alabama**: Constrained optimization (for MBCO extensions)

**Suggested:**

- **MASS**: `mvrnorm()` for parametric bootstrap
- **lavaan**: SEM model support
- **mediation**: mediation package integration
- **medflex**: Natural effects framework support

### Explicit Namespacing

**CRITICAL**: All non-base functions MUST use explicit namespacing:

``` r
# CORRECT
stats::glm(formula, data = data)
stats::coef(model)
lavaan::sem(model, data = data)
MASS::mvrnorm(n, mu, Sigma)

# INCORRECT
glm(formula, data = data)
coef(model)
sem(model, data = data)
mvrnorm(n, mu, Sigma)
```

This prevents namespace conflicts and makes dependencies explicit.

## Important Implementation Details

### P_med Definition

For treatment X, mediator M, and outcome Y:

$$P_{\text{med}} = P\left( Y_{x^{*},M_{x}} > Y_{x,M_{x}} \right)$$

where:

- $Y_{x^{*},M_{x}}$: Outcome under control ($x^{*}$) with mediator at
  its treated level ($M_{x}$)
- $Y_{x,M_{x}}$: Outcome under treatment ($x$) with mediator at its
  treated level ($M_{x}$)

**Interpretation**: Probability that setting treatment to control (while
keeping mediator at treated level) would increase the outcome. Values \>
0.5 indicate positive mediation.

### P_med Computation Methods

**1. Plugin Estimator** (`method = "plugin"`):

- Point estimate only, no confidence interval
- Fast computation via Monte Carlo integration
- Use for quick checks or when inference not needed

**2. Parametric Bootstrap** (`method = "parametric_bootstrap"`):

- Samples from multivariate normal of parameter estimates
- Faster than nonparametric; requires normality assumption
- Default for most applications
- Use: `n_boot = 1000` (minimum), `n_boot = 5000` (recommended)

**3. Nonparametric Bootstrap** (`method = "nonparametric_bootstrap"`):

- Resamples data and refits models
- Robust to parameter distribution assumptions
- Computationally intensive
- Use when sample size adequate (n \> 100) and parametric assumptions
  questionable

### Bootstrap Implementation

Both bootstrap methods follow this pattern:

1.  **Generate bootstrap samples**
    - Parametric: Sample from
      $N\left( \widehat{\theta},\widehat{\Sigma} \right)$
    - Nonparametric: Resample data with replacement
2.  **Compute P_med for each sample**
    - Parametric: Use sampled parameters directly
    - Nonparametric: Refit models, extract parameters, compute P_med
3.  **Extract empirical quantiles**
    - CI: quantiles at $\alpha/2$ and $1 - \alpha/2$
    - Point estimate: mean or median of bootstrap distribution
4.  **Return results**
    - Full bootstrap distribution stored in `@boot_estimates`
    - Enables diagnostic plots and summary statistics

**Memory considerations:**

- Bootstrap distributions stored in memory
- Large `n_boot` with many parameters can consume substantial memory
- For `n_boot > 10000`, consider saving to disk and streaming

### Handling Interactions

When exposure-mediator interactions are present ($X \times M$ in outcome
model):

**Conditioning approach:**

``` r
pmed(
  Y ~ X + M + X:M + C,
  formula_m = M ~ X + C,
  data = data,
  treatment = "X",
  mediator = "M",
  condition_on_x = 0  # Evaluate P_med at X = 0
)
```

**Options for `condition_on_x`:**

- Numeric value: Condition on specific X level
- `"mean"`: Marginalize over X distribution
- Vector: Multiple conditioning values (returns list)

**Mathematical note**: With interactions, P_med becomes conditional on X
level, as the mediation effect varies by treatment level.

### Multiple Mediators

**Parallel mediators** (M1 and M2 both mediate X → Y):

``` r
pmed(
  formula_y = Y ~ X + M1 + M2 + C,
  formula_m = list(
    M1 ~ X + C,
    M2 ~ X + C
  ),
  data = data,
  treatment = "X",
  mediator = c("M1", "M2"),
  effect = "total"  # Total indirect effect
)
```

**Sequential mediators** (X → M1 → M2 → Y):

- Requires path-specific effect definition
- Currently in development (Phase 3)
- Will use nested counterfactuals

## Statistical Assumptions & Diagnostics

### Key Assumptions

1.  **Correct model specification**:
    - Mediator model: E\[M\|X,C\] correctly specified
    - Outcome model: E\[Y\|X,M,C\] correctly specified
    - Misspecification biases P_med estimate
2.  **No unmeasured confounding** (Imai et al., 2010):
    - $\{ Y_{x,m},M_{x}\}\bot X|C$ (exposure-outcome)
    - $Y_{x,m}\bot M|X,C$ (mediator-outcome)
    - $Y_{x,m}\bot M_{x^{*}}$ (cross-world independence)
    - These are *untestable* assumptions; require subject-matter
      knowledge
3.  **Parameter normality** (for parametric bootstrap):
    - $\left( \widehat{a},\widehat{b},\widehat{c}\prime \right) \sim N(\theta,\Sigma)$
    - Generally reasonable for large samples (CLT)
    - Check with Q-Q plots of bootstrap samples
4.  **Correct link functions**:
    - GLM families must match data types
    - Binary outcome: use `family = binomial()`
    - Count outcome: use `family = poisson()`

### Diagnostics

**Check bootstrap distribution:**

``` r
result <- pmed(..., method = "parametric_bootstrap", n_boot = 1000)

# Plot bootstrap distribution
plot(result)

# Summary statistics
summary(result)

# Extract bootstrap samples for custom diagnostics
boot_dist <- result@boot_estimates
hist(boot_dist)
```

**Warning signs:**

- P_med outside \[0, 1\]: Model misspecification likely
- Highly skewed bootstrap distribution: Consider nonparametric bootstrap
- Wide CI relative to estimate: Large uncertainty, increase sample size
- CI includes 0.5: Cannot conclude mediation direction

### Common Pitfalls

1.  **Ignoring confounders**: Always adjust for common causes of M-Y
2.  **Assuming causality from association**: P_med quantifies
    *association*, not *causation*
3.  **Misinterpreting P_med = 0.5**: This is the null (no mediation),
    not “50% mediation”
4.  **Using plugin estimator for inference**: Plugin gives no CI; always
    bootstrap
5.  **Insufficient bootstrap samples**: Use at least n_boot = 1000,
    preferably 5000+

## S7 Method Dispatch

### Adding New Extraction Methods

To add support for a new package (e.g., `newpackage`):

1.  **Define S7 method**:

``` r
#' Extract from newpackage
#'
#' @export
S7::method(extract_mediation, newpackage::fitted_model) <- function(object,
                                                                     treatment,
                                                                     mediator,
                                                                     data = NULL,
                                                                     ...) {
  # Extract parameters
  params <- newpackage::get_params(object)
  
  # Extract a, b, c' paths
  a_path <- params[treatment, mediator]
  b_path <- params[mediator, outcome]
  c_prime <- params[treatment, outcome]
  
  # Get variance-covariance matrix
  vcov_mat <- newpackage::get_vcov(object)
  
  # Return MediationExtract
  MediationExtract(
    estimates = params,
    a_path = a_path,
    b_path = b_path,
    c_prime = c_prime,
    vcov = vcov_mat,
    data = data,
    n_obs = newpackage::nobs(object),
    source_package = "newpackage",
    converged = newpackage::converged(object),
    treatment = treatment,
    mediator = mediator,
    outcome = newpackage::get_outcome(object)
  )
}
```

2.  **Add to DESCRIPTION**:

&nbsp;

    Suggests:
        newpackage

3.  **Add conditional requirement check**:

``` r
if (!requireNamespace("newpackage", quietly = TRUE)) {
  stop("Package 'newpackage' needed. Install with: install.packages('newpackage')")
}
```

### Creating New S7 Classes

When extending functionality with new classes:

``` r
#' New Class for Extension
#'
#' @export
NewClass <- S7::new_class(
  "NewClass",
  package = "probmed",
  parent = MediationExtract,  # Inherit from base class
  
  properties = list(
    # New properties specific to this class
    new_property = S7::class_numeric
  ),
  
  validator = function(self) {
    # Validation logic
    if (self@new_property < 0) {
      "new_property must be non-negative"
    }
  }
)
```

## Testing Strategy

### Unit Tests Should Cover

1.  **P_med Computation Accuracy**
    - Compare to hand-calculated examples
    - Verify P_med ∈ \[0, 1\]
    - Check symmetry properties
    - Test with known effect sizes
2.  **Bootstrap Methods**
    - Parametric vs nonparametric convergence with large n_boot
    - Reproducibility with set.seed()
    - CI coverage in simulations (should be ~95% for 95% CI)
    - Memory efficiency with large n_boot
3.  **S7 Validation**
    - Property type checking works
    - Validators catch invalid inputs
    - Class inheritance functions correctly
4.  **Model Extraction**
    - Identical results from lavaan object vs manual specification
    - Correct extraction from different model types
    - Proper handling of constrained parameters
    - Error handling for convergence failures
5.  **Edge Cases**
    - Small sample sizes (n \< 50)
    - Near-zero mediation effects
    - Perfect mediation (P_med ≈ 1)
    - Suppression effects (P_med \< 0.5)
    - Missing data handling
    - Singular covariance matrices
6.  **Formula Interface**
    - Correct parsing of formulas with interactions
    - Multiple mediators specification
    - Interaction detection and conditioning

### Test Organization

    tests/
    ├── testthat/
    │   ├── test-classes.R           # S7 class validation
    │   ├── test-pmed-glm.R          # GLM interface
    │   ├── test-pmed-bootstrap.R    # Bootstrap methods
    │   ├── test-extract-lavaan.R    # lavaan extraction
    │   ├── test-interactions.R      # Interaction models
    │   ├── test-multiple.R          # Multiple mediators
    │   └── test-edge-cases.R        # Edge cases and errors
    └── manual/
        └── validation-simulations.R  # Long-running validation

### Coverage Expectations

- **Target**: \>90% test coverage
- **Critical paths**: 100% coverage for P_med computation core
- **Bootstrap**: Coverage assessment via simulation studies
- **CI coverage**: Should be within ±2% of nominal level in simulations

## Key Mathematical Formulas

### P_med Definition

$$P_{\text{med}} = P\left( Y_{x^{*},M_{x}} > Y_{x,M_{x}} \right) = \int P\left( Y_{x^{*},M_{x}} > Y_{x,M_{x}}|C = c \right)dF_{C}(c)$$

### Simple Mediation (Linear Models)

For linear outcome model
$Y = \beta_{0} + \beta_{1}X + \beta_{2}M + \beta_{3}C + \epsilon_{Y}$:

$$P_{\text{med}} = \Phi\left( \frac{ab\left( x - x^{*} \right)}{\sqrt{\sigma_{Y}^{2}}} \right)$$

where:

- $a$ = effect of X on M
- $b$ = effect of M on Y
- $\Phi$ = standard normal CDF
- $\sigma_{Y}^{2}$ = residual variance of Y

### Approximate Standard Error (Delta Method)

For parametric bootstrap, standard error approximated via:

$$SE\left( P_{\text{med}} \right) \approx \sqrt{\text{Var}_{\widehat{\theta}}\left\lbrack P_{\text{med}}\left( \widehat{\theta} \right) \right\rbrack}$$

where variance computed empirically from bootstrap samples.

### With Interactions

When
$Y = \beta_{0} + \beta_{1}X + \beta_{2}M + \beta_{3}XM + \beta_{4}C + \epsilon$:

$$P_{\text{med}}\left( x_{0} \right) = P\left( Y_{x^{*},M_{x}} > Y_{x,M_{x}}|X = x_{0} \right)$$

Conditioning value $x_{0}$ must be specified; P_med varies by treatment
level.

## Integration with medrobust Package

The probmed package is designed to work seamlessly with medrobust for
sensitivity analysis:

``` r
library(probmed)
library(medrobust)

# Compute P_med
result <- pmed(
  Y ~ X + M + C,
  formula_m = M ~ X + C,
  data = data,
  treatment = "X",
  mediator = "M"
)

# Sensitivity analysis for unmeasured confounding
sens <- medrobust::sensitivity_analysis(
  result,
  rho_range = seq(-0.5, 0.5, by = 0.1)
)

plot(sens)
```

**medrobust should provide:**

- S7 method for `PmedResult` class
- Sensitivity parameters: ρ (U → M, U → Y correlation)
- Tipping point analysis
- Visualization of P_med under confounding scenarios

## Additional Resources

### Key Publications

**P_med Method:**

- Tofighi et al. (2025). P_med: A scale-free effect size for mediation
  analysis. *Psychological Methods*. (In preparation)

**Causal Mediation Framework:**

- Imai, K., Keele, L., & Tingley, D. (2010). A general approach to
  causal mediation analysis. *Psychological Methods*, 15(4), 309–334.
- VanderWeele, T. J. (2015). *Explanation in Causal Inference: Methods
  for Mediation and Interaction*. Oxford University Press.
- Pearl, J. (2001). Direct and indirect effects. *Proceedings of UAI*,
  411–420.

**Scale-Free Effect Sizes:**

- McGraw, K. O., & Wong, S. P. (1992). A common language effect size
  statistic. *Psychological Bulletin*, 111(2), 361–365.
- Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the
  CL common language effect size statistics of McGraw and Wong. *Journal
  of Educational and Behavioral Statistics*, 25(2), 101–132.

**Multiple Mediators:**

- VanderWeele, T. J., & Vansteelandt, S. (2014). Mediation analysis with
  multiple mediators. *European Journal of Epidemiology*, 29(11),
  801–810.
- Daniel, R. M., De Stavola, B. L., Cousens, S. N., & Vansteelandt, S.
  (2015). Causal mediation analysis with multiple mediators.
  *Biometrics*, 71(1), 1–14.

**Sensitivity Analysis:**

- Imai, K., Keele, L., & Yamamoto, T. (2010). Identification, inference
  and sensitivity analysis for causal mediation effects. *Statistical
  Science*, 25(1), 51–71.
- VanderWeele, T. J. (2010). Bias formulas for sensitivity analysis for
  direct and indirect effects. *Epidemiology*, 21(4), 540–551.

### Package Website

- GitHub: <https://github.com/data-wise/probmed>
- Documentation: <https://data-wise.github.io/probmed>
- CRAN: (Upon release)

### Related Packages

- **mediation**: General mediation analysis (Tingley et al.)
- **lavaan**: Structural equation modeling (Rosseel)
- **medflex**: Natural effects models (Steen et al.)
- **medrobust**: Sensitivity analysis (Tofighi, companion package)
- **RMediation**: Confidence intervals for indirect effects (Tofighi &
  MacKinnon)

## Development Roadmap

### Phase 1: Core Functionality (Current)

S7 class architecture

Formula interface with GLM

Bootstrap inference methods

Basic extraction from lm/glm

Print/summary/plot methods

### Phase 2: Model Integration (Weeks 6-9)

lavaan extraction

mediation package extraction

medflex extraction

Comprehensive testing

### Phase 3: Extensions (Months 3-6)

Exposure-mediator interactions

Multiple parallel mediators

Sequential mediators (path-specific)

Integration with medrobust

### Phase 4: Advanced Features (Future)

Longitudinal mediation

High-dimensional mediators

Bayesian estimation (via Stan/JAGS)

Multiply robust estimation
