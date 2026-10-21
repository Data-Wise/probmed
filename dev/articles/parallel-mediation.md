# Parallel mediation: joint P_med

## Overview

When a treatment influences an outcome through **several mediators at
once** — each modelled in its own equation but all entering a single
outcome equation — the structure is *parallel mediation*:

``` math
X \rightarrow M_1, \quad X \rightarrow M_2, \quad \dots, \quad X \rightarrow M_k,
\qquad Y \sim X + M_1 + \dots + M_k .
```

`probmed` extends $`P_{med}`$ to this setting with a **joint** estimand:
the probability that the outcome is higher when *all* mediators are set
to their treated levels than when all are set to their control levels.
This vignette derives the estimand, gives its closed form, and
demonstrates the four inference methods. It covers the **Gaussian** case
(continuous mediators and outcome).

## The joint estimand

Write the treated and control treatment levels as \$x = \$ `x_value` and
\$x^= \$ `x_ref`, with contrast $`\delta = x - x^\ast`$. As in the
single-mediator case, the treatment is **held at $`x`$ for both**
potential outcomes; only the mediators vary between their treated levels
$`M_j(x)`$ and control levels $`M_j(x^\ast)`$, drawn independently. The
direct effect $`c'`$ therefore cancels and $`P_{med}`$ reflects
mediation only.

Each mediator follows $`M_j = i_{m_j} + a_j X + \varepsilon_{m_j}`$ with
$`\varepsilon_{m_j} \sim N(0, \sigma_{m_j}^2)`$, and the outcome follows
$`Y = i_y + c' X + \sum_j b_j M_j + \varepsilon_y`$ with
$`\varepsilon_y \sim N(0, \sigma_Y^2)`$. The two potential outcomes
(treatment fixed at $`x`$) are

``` math
Y_t = i_y + c' x + \sum_j b_j M_j(x) + \varepsilon_t,
\qquad
Y_c = i_y + c' x + \sum_j b_j M_j(x^\ast) + \varepsilon_c .
```

Each mediator difference is Gaussian and independent across $`j`$,

``` math
M_j(x) - M_j(x^\ast) \sim N\!\big(a_j\delta,\; 2\sigma_{m_j}^2\big),
```

so the outcome difference is Gaussian with

``` math
Y_t - Y_c \;\sim\; N\!\Big(\; \delta \textstyle\sum_j a_j b_j,\;\;
2\textstyle\sum_j b_j^2 \sigma_{m_j}^2 + 2\sigma_Y^2 \;\Big).
```

The **joint parallel $`P_{med}`$** is then a clean closed form:

``` math
\boxed{\;P_{med} = \Phi\!\left(
\frac{\delta \sum_j a_j b_j}
     {\sqrt{2\sum_j b_j^2 \sigma_{m_j}^2 + 2\sigma_Y^2}}
\right)\;}
```

which recovers the single-mediator formula
$`\Phi\!\big(ab / \sqrt{2(b^2\sigma_M^2 + \sigma_Y^2)}\big)`$ exactly at
$`k = 1`$. Reversing the contrast ($`\delta \to -\delta`$) maps
$`P_{med} \to 1 - P_{med}`$. The associated **total indirect effect** is
$`\sum_j a_j b_j`$.

## A worked example

We simulate two parallel mediators, fit each equation, and extract a
`ParallelMediationData` object with **medfit**.

``` r

n <- 1500
X <- rnorm(n)
M1 <- 0.5 * X + rnorm(n)
M2 <- 0.4 * X + rnorm(n)
Y <- 0.3 * X + 0.5 * M1 + 0.6 * M2 + rnorm(n)
dat <- data.frame(X, M1, M2, Y)

fit_m1 <- lm(M1 ~ X, data = dat)
fit_m2 <- lm(M2 ~ X, data = dat)
fit_y <- lm(Y ~ X + M1 + M2, data = dat)

ex <- medfit::extract_mediation(
  object = fit_m1, model_y = fit_y,
  treatment = "X", mediator = c("M1", "M2"),
  mediator_models = list(fit_m2),
  structure = "parallel", data = dat
)
```

### Plugin (point estimate)

``` r

library(probmed)
pmed(ex, method = "plugin")
#> 
#> P_med: Probability of Mediated Shift
#> ====================================
#> 
#> Estimate: 0.608 
#> 
#> Indirect Effect (Product of Coefficients):
#> Estimate: 0.508 
#> 
#> Inference: plugin 
#> 
#> Treatment contrast: X = 1 vs. X* = 0 
#> 
#> Interpretation:
#>   P(Y(1, M(1)) > Y(1, M(0))) = 0.608
#>   P that the mediator shift (M(0) -> M(1)) leaves a random individual better off, holding X = 1.
```

The point estimate is the closed form above; the indirect effect is
$`\sum_j a_j b_j`$.

### Parametric bootstrap

Draws the structural coefficients $`(a_1, \dots, b_k)`$ from
$`N(\hat\theta, \widehat{\mathrm{Var}})`$ (residual SDs held fixed) and
evaluates the closed form per draw.

``` r

pmed(ex, method = "parametric_bootstrap", n_boot = 500, seed = 1)
#> 
#> P_med: Probability of Mediated Shift
#> ====================================
#> 
#> Estimate: 0.608 
#> 95% CI: [ 0.599 ,  0.618 ]
#> 
#> Indirect Effect (Product of Coefficients):
#> Estimate: 0.508 
#> 95% CI: [ 0.455 ,  0.563 ]
#> 
#> Inference: parametric_bootstrap 
#> Bootstrap samples: 500 
#> 
#> Treatment contrast: X = 1 vs. X* = 0 
#> 
#> Interpretation:
#>   P(Y(1, M(1)) > Y(1, M(0))) = 0.608
#>   P that the mediator shift (M(0) -> M(1)) leaves a random individual better off, holding X = 1.
```

### Nonparametric bootstrap

Resamples rows, refits the mediator and outcome models, and evaluates
the closed form per resample.

``` r

pmed(ex, method = "nonparametric_bootstrap", n_boot = 300, seed = 1)
#> 
#> P_med: Probability of Mediated Shift
#> ====================================
#> 
#> Estimate: 0.608 
#> 95% CI: [ 0.599 ,  0.618 ]
#> 
#> Indirect Effect (Product of Coefficients):
#> Estimate: 0.506 
#> 95% CI: [ 0.457 ,  0.562 ]
#> 
#> Inference: nonparametric_bootstrap 
#> Bootstrap samples: 300 
#> 
#> Treatment contrast: X = 1 vs. X* = 0 
#> 
#> Interpretation:
#>   P(Y(1, M(1)) > Y(1, M(0))) = 0.608
#>   P that the mediator shift (M(0) -> M(1)) leaves a random individual better off, holding X = 1.
```

### MBCO (deterministic)

Inverts a constrained likelihood-ratio test for both the joint
$`P_{med}`$ and the total indirect effect — no resampling, so it is
seed-free and grid-resolution-independent.

``` r

pmed(ex, method = "mbco")
#> 
#> P_med: Probability of Mediated Shift
#> ====================================
#> 
#> Estimate: 0.608 
#> 95% CI: [ 0.599 ,  0.619 ]
#> 
#> Indirect Effect (Product of Coefficients):
#> Estimate: 0.508 
#> 95% CI: [ 0.456 ,  0.563 ]
#> 
#> Inference: mbco 
#> 
#> Treatment contrast: X = 1 vs. X* = 0 
#> 
#> Interpretation:
#>   P(Y(1, M(1)) > Y(1, M(0))) = 0.608
#>   P that the mediator shift (M(0) -> M(1)) leaves a random individual better off, holding X = 1.
```

## Scope and limitations

- **Gaussian only.** Parallel $`P_{med}`$ requires continuous mediators
  and a continuous outcome. A non-Gaussian `ParallelMediationData` (one
  lacking finite residual SDs) raises an informative error; non-Gaussian
  support is blocked pending family slots in the medfit parallel
  container.
- **Joint estimand only.** This release reports the single joint
  $`P_{med}`$ for moving all mediators together. Path-specific
  ($`P_{med}`$ per mediator) decomposition is planned future work.
- **Object entry.** Reach the method via
  `medfit::extract_mediation(..., structure = "parallel")` then
  [`pmed()`](https://data-wise.github.io/probmed/dev/reference/pmed.md);
  a formula convenience wrapper is planned.
- **Serial / sequential** mediation ($`X \to M_1 \to M_2 \to Y`$) is a
  separate estimand and is not covered here.
