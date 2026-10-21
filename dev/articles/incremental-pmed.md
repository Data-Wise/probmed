# Incremental mediated elasticity: the P_med^delta curve

## Why a curve, not a number

The standard proportion mediated $`P_{med} = \text{IIE} / \text{OE}`$
answers the question “what share of the effect passes through the
mediator?” for one fixed contrast $`A = 0`$ vs $`A = 1`$. But that
contrast is arbitrary when treatment is a continuous or graded dose. The
**incremental mediated elasticity** replaces the global contrast with a
derivative question: if assignment odds are tilted by a factor
$`\delta`$, what fraction of the resulting marginal shift in the outcome
travels through the mediator?

Formally, the **tilted assignment probability** is

``` math
q_\delta(C) = \frac{\delta \, g(C)}{\delta \, g(C) + 1 - g(C)},
\qquad g(C) = P(A = 1 \mid C),
```

and the estimand is

``` math
P_{med}^\delta(\delta) = \frac{\theta'_{med}(\delta)}{\theta'_{tot}(\delta)},
```

where $`\theta'_{dir}`$ and $`\theta'_{med}`$ are the marginal
(derivative-scale) direct and mediated elasticities, and
$`\theta'_{tot} = \theta'_{dir} + \theta'_{med}`$ is their exact sum (no
interaction remainder $`R(\delta) = 0`$ for all $`\delta`$). The curve
$`P_{med}^\delta(\delta)`$ is the quantity that varies with $`\delta`$:
it is **flat** under a linear-Gaussian model with no $`A \times M`$
interaction, and it **bends** — revealing the interaction — when one is
present.

## A worked example

The DGP below places an $`A \times M`$ interaction in the outcome so the
curve is non-flat.
[`incr_pmed()`](https://data-wise.github.io/probmed/dev/reference/incr_pmed.md)
evaluates the estimand at five tilt values and returns an
`IncrPmedResult` whose print method displays one row per $`\delta`$.

``` r

library(probmed)
```

``` r

n <- 1500
C <- rnorm(n)
A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
M <- 0.6 * A + 0.4 * C + rnorm(n)
Y <- 0.5 * A + 0.7 * M + 0.5 * A * M + 0.3 * C + rnorm(n)   # A*M interaction
d <- data.frame(A, M, Y, C)
fit <- incr_pmed(d, deltas = c(1/3, 1/2, 1, 2, 3))
fit
#> Incremental mediated elasticity P_med^delta(delta) (onestep-crossfit, n=1500)
#>   delta=0.333  P_med=0.490  [0.404, 0.577]  (dir=0.288 med=0.276 tot=0.564)
#>   delta=0.5    P_med=0.494  [0.415, 0.574]  (dir=0.223 med=0.218 tot=0.441)
#>   delta=1      P_med=0.501  [0.428, 0.573]  (dir=0.129 med=0.130 tot=0.259)
#>   delta=2      P_med=0.505  [0.436, 0.575]  (dir=0.065 med=0.067 tot=0.132)
#>   delta=3      P_med=0.507  [0.438, 0.577]  (dir=0.041 med=0.042 tot=0.083)
```

Each row of the printed table gives the mediated elasticity (`med`),
total elasticity (`tot`), and $`P_{med}^\delta`$ with its confidence
interval at the corresponding $`\delta`$. Because the interaction term
amplifies the mediated share as $`\delta`$ grows (larger tilts assign
more treated units and the $`A \times M`$ term compounds), the curve
rises rather than remaining flat. A flat curve — constant `Pmed` across
$`\delta`$ — certifies that no interaction is present and the classical
proportion mediated tells the whole story. A bending curve tells you
that the “share through the mediator” depends on how aggressively you
shift treatment, which a single number would obscure.

## Inference

The estimator behind each row is a **cross-fitted one-step** (Kennedy
2019). The propensity score $`g(C)`$ is fitted in cross-fitting folds
and the pseudo-outcomes are constructed to be Neyman-orthogonal in
$`g`$: the influence function includes the full g-score (Kennedy 2019,
Corollary 2, Term II). Orthogonality in the propensity means that
first-order errors in $`\hat g`$ do not inflate the asymptotic variance
of the elasticity estimates — inference remains valid at standard
nonparametric rates even when $`g`$ is estimated by a flexible learner.

Note carefully what orthogonality does *not* imply: because
$`P_{med}^\delta`$ is itself a functional of $`g`$ (the tilt
$`q_\delta`$ is a transformation of $`g`$), **double robustness of the
point estimate is impossible** (Kennedy 2019, §3.3). Misspecifying $`g`$
biases $`\hat q_\delta`$, which biases the elasticity estimates. The
efficiency and inferential robustness that orthogonality buys are
genuine; DR of the point estimate is not on offer for this estimand, and
probmed does not claim it.

Pointwise confidence intervals at each $`\delta`$ are computed via the
ratio-identity influence function for $`\text{med} / \text{tot}`$.
Intervals are printed as `[lo, hi]` at the `ci_level` stored in the
result object (default 95 %).

## Sensitivity to M–Y confounding

Unmeasured confounding between the mediator and outcome is the canonical
threat to mediation conclusions.
[`incr_sensitivity()`](https://data-wise.github.io/probmed/dev/reference/incr_sensitivity.md)
applies the additive-offset sensitivity model of
[`pmed_sensitivity()`](https://data-wise.github.io/probmed/dev/reference/pmed_sensitivity.md)
to every row of the curve: it reports, for each $`\delta`$, the M–Y bias
magnitude that would drive $`P_{med}^\delta`$ to zero.

``` r

incr_sensitivity(fit, threshold = 0)
#>       delta        med        tot      Pmed     tipping tipping_threshold
#> 1 0.3333333 0.27643023 0.56394258 0.4901744 -0.27643023       -0.27643023
#> 2 0.5000000 0.21822093 0.44135523 0.4944338 -0.21822093       -0.21822093
#> 3 1.0000000 0.12968367 0.25906367 0.5005861 -0.12968367       -0.12968367
#> 4 2.0000000 0.06667354 0.13193290 0.5053594 -0.06667354       -0.06667354
#> 5 3.0000000 0.04224607 0.08331181 0.5070838 -0.04224607       -0.04224607
```

The column `tipping` is $`-\theta'_{med}(\delta)`$ — the additive
numerator bias that would make the mediated share vanish at this
$`\delta`$. `tipping_threshold` gives the corresponding quantity for the
user-supplied `threshold` (here 0, so both columns coincide). Larger
`|tipping|` at a given $`\delta`$ signals that the mediated share at
that tilt is more durable to unmeasured M–Y confounding. When the curve
rises (as here, due to the interaction), the rows at large $`\delta`$
typically have the largest absolute tipping values, meaning the
interaction-driven amplification is also the hardest to explain away.

## Summary

- [`incr_pmed()`](https://data-wise.github.io/probmed/dev/reference/incr_pmed.md)
  estimates the incremental mediated elasticity curve
  $`P_{med}^\delta(\delta)`$ at a user-supplied grid of tilt factors.
- A flat curve signals no $`A \times M`$ interaction; a bending curve
  reveals one.
- Inference is Neyman-orthogonal in the propensity but not doubly-robust
  in the point estimate (the estimand is itself a functional of $`g`$).
- [`incr_sensitivity()`](https://data-wise.github.io/probmed/dev/reference/incr_sensitivity.md)
  extends the curve object with per-$`\delta`$ tipping-point analyses
  for M–Y unmeasured confounding.
