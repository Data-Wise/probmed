# Incremental Mediated Elasticity Curve

Estimate the incremental mediated elasticity `P_med^delta(delta)` – the
derivative-scale proportion mediated as a function of the treatment-tilt
factor `delta`. The estimator is cross-fitted and built on the same
triply-robust corner pseudo-outcomes as
[`ward_residual()`](https://data-wise.github.io/probmed/reference/ward_residual.md);
standard errors use the ratio-identity influence function for
`med / tot` at each `delta`.

## Usage

``` r
incr_pmed(
  object,
  deltas = c(1/3, 1/2, 1, 2, 3),
  covars = "C",
  K = 5L,
  ci_level = 0.95,
  seed = 1L,
  ...
)
```

## Arguments

- object:

  A data frame with columns `A` (binary treatment), `M` (mediator), `Y`
  (outcome), and the covariates named in `covars`.

- deltas:

  Numeric: tilt factors at which to evaluate the curve.

- covars:

  Character: covariate column names (default `"C"`).

- K:

  Integer: number of cross-fitting folds (default `5`).

- ci_level:

  Numeric: confidence level (default `0.95`).

- seed:

  Integer: random seed for fold assignment (default `1`).

- ...:

  Unused.

## Value

An
[IncrPmedResult](https://data-wise.github.io/probmed/reference/IncrPmedResult.md)
object whose `@curve` slot holds one row per `delta` with the direct,
mediated, and total elasticities, the proportion mediated `Pmed`, and
its standard error and Wald interval.

## Details

At each `delta` the tilted assignment probability is
`q = delta * g / (delta * g + 1 - g)` with derivative
`q' = g (1 - g) / (delta * g + 1 - g)^2`, where `g = P(A = 1 | C)`. The
direct and mediated elasticities are weighted contrasts of the four
corner pseudo-outcomes by `q` and `q'`; by Theorem 1 (multivariate chain
rule) they sum to the total elasticity exactly, so
`P_med^delta = med / (dir + med)`.

**Inference.** Point estimates use the tilt-derivative weight
`q' = g(1 - g) / (delta g + 1 - g)^2` (Kennedy 2019, Corollary 2, Term
I). Standard errors use the full efficient influence function for the
ratio `P_med = med / tot`, which adds the g-score correction
`(dq/dg) * (A - g(C)) * (dir * E[a_med | C] - med * E[a_dir | C]) / tot^2`,
where `a_med = q' * gamma_med` and `a_dir = q' * gamma_dir` are the same
tilt-weighted corner contrasts used in the point estimate (the g-score
term must carry the `q'` weight to match the units of `dir`/`med`/`tot`)
and `E[. | C]` is a linear projection on the covariates, with
`dq/dg = delta / (delta g + 1 - g)^2` (Term II). Term II is mean-zero
(by `E[A - g(C) | C] = 0`), so it does not shift the point estimate but
restores Neyman orthogonality w.r.t. the propensity score, ensuring the
CI is consistent under nonparametric estimation of `g`. The projection
matters: the per-row `a_med`/`a_dir` carry the corner EIF's
inverse-probability noise, and multiplying that by `A - g(C)` inflated
the se by 1.40x at `delta = 2` (fixed 2026-08-22; se/empSD now 0.94-0.97
across `delta`).

**Scope (exposure type).** The incremental intervention is defined only
for a **binary** treatment `A`: the odds-ratio tilt
`q = delta g / (delta g + 1 - g)` requires `g = P(A = 1 | C)` to be a
Bernoulli probability. Continuous or multi-valued exposures are a
distinct object requiring a *modified treatment policy* (shift)
intervention with its own estimand and efficient influence function;
they are deliberately out of scope here (deferred to future work – a
separate MTP-mediation paper).

## Examples

``` r
set.seed(1)
n <- 1500
C <- rnorm(n)
A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
M <- 0.6 * A + 0.4 * C + rnorm(n)
Y <- 0.5 * A + 0.7 * M + 0.3 * C + rnorm(n)
d <- data.frame(A, M, Y, C)
incr_pmed(d, deltas = c(0.5, 1, 2))
#> Incremental mediated elasticity P_med^delta(delta) (onestep-crossfit, n=1500)
#>   delta=0.5    P_med=0.499  [0.405, 0.593]  (dir=0.164 med=0.163 tot=0.327)
#>   delta=1      P_med=0.493  [0.405, 0.580]  (dir=0.097 med=0.094 tot=0.191)
#>   delta=2      P_med=0.486  [0.401, 0.570]  (dir=0.049 med=0.046 tot=0.095)
```
