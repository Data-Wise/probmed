# Introduction to probmed

## Overview

The `probmed` package provides a framework for computing
**$`P_{med}`$**, a scale-free probabilistic effect size for causal
mediation analysis. Unlike traditional effect sizes like the “Proportion
Mediated” ($`P_M`$), which can be unstable or difficult to interpret
(especially with inconsistent mediation or non-linear models),
$`P_{med}`$ offers a clear probabilistic interpretation.

### What is $`P_{med}`$?

$`P_{med}`$ is defined as the probability that the outcome for an
individual in the treatment group ($`Y(1)`$) is superior to the outcome
for an individual in the control group ($`Y(0)`$), specifically due to
the indirect effect transmitted through the mediator ($`M`$).

Mathematically, for a treatment $`X`$, mediator $`M`$, and outcome
$`Y`$:

``` math
 P_{med} = P(Y(1, M(1)) > Y(0, M(0))) 
```

where the direct effect of $`X`$ on $`Y`$ is held constant or accounted
for, focusing purely on the mediation path.

Key advantages:

1.  **Scale-Free**: It does not depend on the units of $`Y`$ or $`M`$.
2.  **Bounded**: It ranges from 0 to 1 (or 0.5 for no effect in some
    contexts).
3.  **Robust**: It works well for both linear and generalized linear
    models (GLMs).

## Installation

You can install the development version of `probmed` from GitHub:

``` r

# install.packages("devtools")
devtools::install_github("data-wise/probmed")
```

## Basic Example: Linear Mediation

Let’s demonstrate `probmed` with a simple linear mediation model. We
will simulate data where $`X`$ affects $`M`$, and both $`X`$ and $`M`$
affect $`Y`$.

### 1. Simulate Data

``` r

library(probmed)

set.seed(123)
n <- 500

# Confounder / Covariate
C <- rnorm(n)

# Treatment (Continuous)
X <- rnorm(n)

# Mediator: M ~ X + C
# Path a = 0.5
M <- 0.5 * X + 0.3 * C + rnorm(n)

# Outcome: Y ~ M + X + C
# Path b = 0.4, Path c' = 0.2
Y <- 0.4 * M + 0.2 * X + 0.2 * C + rnorm(n)

data <- data.frame(X, M, Y, C)
head(data)
```

                X          M             Y           C
    1 -0.60189285 -1.4648878 -1.6394155324 -0.56047565
    2 -0.99369859 -1.6058576 -1.1943754836 -0.23017749
    3  1.02678506  0.9630248  0.0002105782  1.55870831
    4  0.75106130  0.2645080  0.8971858963  0.07050839
    5 -1.50916654 -3.2651397 -0.4616766212  0.12928774
    6 -0.09514745  1.5075192  3.0542047502  1.71506499

### 2. Estimate $`P_{med}`$

To estimate $`P_{med}`$, we use the
[`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md)
function. We need to specify:

- The **outcome model** formula (`Y ~ ...`).
- The **mediator model** formula (`formula_m = ...`).
- The names of the **treatment** and **mediator** variables.

We will use the **parametric bootstrap** method for inference, which is
fast and accurate for standard models.

``` r

result <- pmed(
    Y ~ X + M + C, # Outcome model
    formula_m = M ~ X + C, # Mediator model
    data = data,
    treatment = "X",
    mediator = "M",
    method = "parametric_bootstrap",
    n_boot = 1000, # Number of bootstrap samples
    seed = 42
)
```

### 3. Interpret Results

We can print and summarize the results.

``` r

print(result)
```


    P_med: Probability of Mediated Shift
    ====================================

    Estimate: 0.35
    95% CI: [ 0.326 ,  0.374 ]

    Inference: parametric_bootstrap
    Bootstrap samples: 1000

    Treatment contrast: X = 1 vs. X* = 0

    Interpretation:
      P(Y_{X*, M_X} > Y_{X, M_X}) = 0.35 

The output shows the estimated $`P_{med}`$ and its 95% confidence
interval.

- **Estimate**: A value significantly greater than 0.50 indicates a
  positive mediation effect.
- **Confidence Interval**: If the interval excludes 0.50, the effect is
  statistically significant.

You can also get a more detailed summary:

``` r

summary(result)
```


    P_med: Probability of Mediated Shift
    ====================================

    Estimate: 0.35
    95% CI: [ 0.326 ,  0.374 ]

    Inference: parametric_bootstrap
    Bootstrap samples: 1000

    Treatment contrast: X = 1 vs. X* = 0

    Interpretation:
      P(Y_{X*, M_X} > Y_{X, M_X}) = 0.35

    Bootstrap Distribution:
       Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
     0.3183  0.3402  0.3494  0.3496  0.3584  0.3904

    Standard Error: 0.01271875

    Source: stats
    Sample size: 500 

## Advanced Example: Binary Outcome (GLM)

`probmed` shines when dealing with non-linear models where traditional
coefficients ($`a \times b`$) are hard to interpret. Let’s look at a
case with a **binary outcome**.

### 1. Simulate Binary Data

``` r

set.seed(456)
n <- 500
X <- rnorm(n)
M <- 0.5 * X + rnorm(n)

# Outcome is binary (Logistic model)
# logit(P(Y=1)) = 0.5*M + 0.2*X
logits <- 0.5 * M + 0.2 * X
probs <- 1 / (1 + exp(-logits))
Y_bin <- rbinom(n, 1, probs)

data_bin <- data.frame(X, M, Y_bin)
```

### 2. Estimate with GLM

We simply specify `family_y = binomial()` to tell `probmed` that the
outcome model is a logistic regression.

``` r

result_bin <- pmed(
    Y_bin ~ X + M,
    formula_m = M ~ X,
    data = data_bin,
    family_y = binomial(), # Specify binomial family
    treatment = "X",
    mediator = "M",
    method = "parametric_bootstrap",
    n_boot = 1000,
    seed = 42
)

print(result_bin)
```


    P_med: Probability of Mediated Shift
    ====================================

    Estimate: 0.378
    95% CI: [ 0.352 ,  0.405 ]

    Inference: parametric_bootstrap
    Bootstrap samples: 1000

    Treatment contrast: X = 1 vs. X* = 0

    Interpretation:
      P(Y_{X*, M_X} > Y_{X, M_X}) = 0.378 

The interpretation remains the same: $`P_{med}`$ is the probability of
the outcome being “better” (higher latent utility or probability class)
due to mediation.

## Conclusion

The `probmed` package makes it easy to compute interpretable, scale-free
effect sizes for mediation analysis. Whether you are working with simple
linear models or complex GLMs, $`P_{med}`$ provides a unified metric for
understanding the magnitude of indirect effects.
