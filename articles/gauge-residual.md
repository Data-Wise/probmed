# Gauge-calibrated proportion mediated: the residual W

## The problem with a single proportion-mediated number

The interventional proportion mediated `P_med = IIE / OE` summarizes
mediation in one number. That summary is only trustworthy when the
overall effect **decomposes additively** into a direct and an indirect
part. With a treatment-by-mediator interaction (or intermediate
confounding) the additive split fails: there is a leftover

``` math
R = \mathrm{OE} - \mathrm{IDE} - \mathrm{IIE},
```

the **non-decomposability remainder**. When $`R`$ is nonzero, `IIE / OE`
is reporting a fraction of a quantity that was never the sum of its
parts.

[`ward_residual()`](https://data-wise.github.io/probmed/reference/ward_residual.md)
estimates `P_med` *and* the **gauge residual**

``` math
W = R / \mathrm{OE},
```

a scale-free flag for that failure. $`W \approx 0`$ certifies the
additive split; a $`W`$ significantly away from 0 says “interpret
`P_med` with care.”

``` r

library(probmed)
```

## A clean (additive) example

With no interaction term, $`R`$ — and hence $`W`$ — should sit at zero,
and the gauge `P_med` agrees with the usual story.

``` r

n <- 1500
C <- rnorm(n)
A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
M <- 0.6 * A + 0.4 * C + rnorm(n)
Y <- 0.5 * A + 0.7 * M + 0.3 * C + rnorm(n)        # no A*M term
ward_residual(data.frame(A, M, Y, C))
#> Gauge-calibrated proportion mediated (onestep-crossfit, n=1500)
#>   P_med = 0.502  Wald [0.310, 0.695]
#>     Fieller 95% CI [0.321, 0.714]
#>   W=R/OE = -0.053  Wald [-0.261, 0.155]  (p=0.616)
#>   OE=0.872  IDE=0.480  IIE=0.438  R=-0.046
```

The `W` confidence interval covers 0: the single-number `P_med` is safe
here.

## An interaction breaks the split

Add an `A * M` term and the gauge lights up — `W` is now significantly
positive and the printout warns that the additive split is unreliable.

``` r

Yi <- 0.5 * A + 0.7 * M + 0.8 * A * M + 0.3 * C + rnorm(n)
ward_residual(data.frame(A, M, Y = Yi, C))
#> Gauge-calibrated proportion mediated (onestep-crossfit, n=1500)
#>   P_med = 0.288  Wald [0.163, 0.413]
#>     Fieller 95% CI [0.166, 0.418]
#>   W=R/OE = 0.427  Wald [0.284, 0.569]  (p=4.78e-09)
#>   OE=1.392  IDE=0.398  IIE=0.401  R=0.594
#>   ! |W| large: additive split unreliable; interpret P_med with care.
```

## Non-binary-coded exposures: `a0` / `a1`

The treatment need not be coded `0/1`. Any **two-level** coding works —
a factor, `{1, 2}`, `{-1, 1}` — by naming the reference `a0` and
comparison `a1`. `A` is recoded internally to the indicator `(A == a1)`.

``` r

dat <- data.frame(A = factor(ifelse(A == 1, "drug", "placebo")),
                  M = M, Y = Yi, C = C)
ward_residual(dat, a0 = "placebo", a1 = "drug")
#> Gauge-calibrated proportion mediated (onestep-crossfit, n=1500)
#>   P_med = 0.288  Wald [0.163, 0.413]
#>     Fieller 95% CI [0.166, 0.418]
#>   W=R/OE = 0.427  Wald [0.284, 0.569]  (p=4.78e-09)
#>   OE=1.392  IDE=0.398  IIE=0.401  R=0.594
#>   ! |W| large: additive split unreliable; interpret P_med with care.
```

`W = R / OE` is **anti-symmetric** in the contrast direction: swapping
`a0` and `a1` flips the sign of `W` (the remainder `R` is invariant,
`OE` changes sign) while leaving its magnitude — the thing you interpret
— unchanged.

A treatment with **more than two** levels is rejected rather than
silently subsetted: restricting the data to two of several levels would
change the population the corner means are averaged over, estimating a
different (sub-population) gauge. Filter to the two levels you intend to
contrast first.

``` r

multi <- data.frame(A = sample(0:2, n, replace = TRUE), M = M, Y = Yi, C = C)
ward_residual(multi, a0 = 0, a1 = 1)
#> Error:
#> ! A has more than two levels (extra: {2}). ward_residual contrasts exactly two levels; subsetting a multi-valued A would change the estimand. Pre-filter to the two levels you intend, or see the manuscript for the multi-valued generalization (not yet implemented).
```

## Inference for a skewed ratio

`W` and `P_med` are **ratios**, with right-skewed sampling
distributions. The default `se_method = "analytic"` reports the
influence-function standard error and a symmetric Wald interval, which
is convenient but mildly anti-conservative for the ratio near the null.
Setting `se_method = "bootstrap"` switches the `W` and `P_med` intervals
to the tail-aware **percentile** bootstrap — the appropriate
construction for a skewed ratio (widening a symmetric SE fails to cover
it). The printed label tracks which construction produced the interval.

``` r

ward_residual(data.frame(A, M, Y = Yi, C), se_method = "bootstrap", B = 200)
#> Gauge-calibrated proportion mediated (onestep-crossfit, n=1500)
#>   P_med = 0.288  percentile [0.213, 0.408]
#>     Fieller 95% CI [0.166, 0.418]
#>   W=R/OE = 0.427  percentile [0.250, 0.510]  (p=5.79e-08)
#>   OE=1.392  IDE=0.398  IIE=0.401  R=0.594
#>   ! |W| large: additive split unreliable; interpret P_med with care.
```

For the near-null / weak-identification regime — `OE` not bounded away
from 0 — prefer the **Fieller** confidence set (`fieller = TRUE`, on by
default), which is honest about being unbounded when the total effect is
not significant, where both Wald and percentile intervals understate the
uncertainty.

## Summary

- [`ward_residual()`](https://data-wise.github.io/probmed/reference/ward_residual.md)
  reports `P_med` together with the gauge residual `W = R/OE`; `W` flags
  when the additive direct/indirect split — and therefore the
  single-number `P_med` — is unreliable.
- Use `a0` / `a1` for any two-level exposure coding; `>2` levels is an
  error by design.
- For ratio inference, `se_method = "bootstrap"` gives percentile
  intervals and the Fieller set handles the near-null case.
