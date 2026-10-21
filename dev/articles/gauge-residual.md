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

[`ward_residual()`](https://data-wise.github.io/probmed/dev/reference/ward_residual.md)
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
#>   P_med = 0.446  Wald [0.355, 0.538]
#>     Fieller 95% CI [0.357, 0.542]
#>   W=R/OE = 0.059  Wald [-0.031, 0.150]  (p=0.199)
#>   OE=0.971  IDE=0.480  IIE=0.433  R=0.058
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
#>   P_med = 0.322  Wald [0.264, 0.380]
#>     Fieller 95% CI [0.264, 0.381]
#>   W=R/OE = 0.444  Wald [0.337, 0.550]  (p=3.41e-16)
#>   OE=1.423  IDE=0.333  IIE=0.459  R=0.632
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
#>   P_med = 0.322  Wald [0.264, 0.380]
#>     Fieller 95% CI [0.264, 0.381]
#>   W=R/OE = 0.444  Wald [0.337, 0.550]  (p=3.41e-16)
#>   OE=1.423  IDE=0.333  IIE=0.459  R=0.632
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

## Inference: what the two arms actually do

`W` and `P_med` are **ratios**. The default `se_method = "analytic"`
reports the influence-function standard error of the cross-fit one-step
estimator and a symmetric Wald interval. Measured against exact truth
(`inst/sim/se_shipped_check.R`, 250 datasets per cell, n = 800) it
covers 0.936 for a continuous outcome and 0.972 for a binary one under
strong identification, 0.912 when the total effect is moderate, and its
dispersion from dataset to dataset matches that of an oracle that knows
the nuisance functions. `reps > 1` averages the influence matrix over
several independent partitions; it removes the fold-split component of
the variance, which with correct weights is small (empSD 0.063 at
`reps = 1` vs 0.061 at `reps = 10`, oracle 0.056):

``` r

ward_residual(data.frame(A, M, Y = Yi, C), reps = 10)
#> Gauge-calibrated proportion mediated (onestep-crossfit, n=1500)
#>   P_med = 0.327  Wald [0.269, 0.384]
#>     Fieller 95% CI [0.269, 0.385]
#>   W=R/OE = 0.437  Wald [0.332, 0.543]  (p=5.12e-16)
#>   OE=1.410  IDE=0.332  IIE=0.461  R=0.617
#>   ! |W| large: additive split unreliable; interpret P_med with care.
```

A note on history, because the earlier documentation said the opposite:
before 2026-08-22 the corner influence function used one
inverse-probability weight per fold —
[`ifelse()`](https://rdrr.io/r/base/ifelse.html) on a scalar test
silently returned the first fold row’s propensity for every row — and
every coverage figure the package reported (Wald ~0.85–0.90, percentile
~1.00, the gate-validation grid, the variance decomposition in
`docs/specs/FINDINGS-2026-08-22-phi-decomposition.md`) was produced with
that bug. Those numbers are superseded; the bootstrap arm and the two
gates were re-measured the same day
(`docs/specs/FINDINGS-2026-08-22-postfix-bootstrap-gates.md`).

Setting `se_method = "bootstrap"` switches the `W` and `P_med` intervals
to **percentile** intervals of a nonparametric bootstrap that resamples
rows and refits the whole cross-fit, with a fresh partition each
resample. With correct weights it buys nothing: percentile coverage
0.905–0.930 in the regular cells at the same width as the Wald interval,
and a bootstrap standard error that is unstable for binary outcomes. It
remains available for users who want a refit bootstrap. The `weak_id`
flag it enables — percentile interval at least 3× wider than Wald — is a
symptom of bootstrap instability near `OE = 0`, not a coverage
diagnostic: flagged and unflagged draws cover alike at every threshold
tried; it no longer warns. The printed label tracks which construction
produced the interval.

``` r

ward_residual(data.frame(A, M, Y = Yi, C), se_method = "bootstrap", B = 200)
#> Gauge-calibrated proportion mediated (onestep-crossfit, n=1500)
#>   P_med = 0.322  percentile [0.266, 0.384]
#>     Fieller 95% CI [0.264, 0.381]
#>   W=R/OE = 0.444  percentile [0.337, 0.535]  (p=2.36e-17)
#>   OE=1.423  IDE=0.333  IIE=0.459  R=0.632
#>   ! |W| large: additive split unreliable; interpret P_med with care.
```

For the near-null / weak-identification regime — `OE` not bounded away
from 0 — the ratio itself is non-regular: in the grid’s near-null cells
even an oracle with known nuisances produces exploding standard errors,
and `oe_regular = FALSE` flags exactly this. Prefer the **Fieller**
confidence set there (`fieller = TRUE`, on by default), which is honest
about being unbounded when the total effect is not significant, and
report the unnormalized effects.

## Summary

- [`ward_residual()`](https://data-wise.github.io/probmed/dev/reference/ward_residual.md)
  reports `P_med` together with the gauge residual `W = R/OE`; `W` flags
  when the additive direct/indirect split — and therefore the
  single-number `P_med` — is unreliable.
- Use `a0` / `a1` for any two-level exposure coding; `>2` levels is an
  error by design.
- For ratio inference, the default analytic Wald interval is calibrated
  under regular identification (0.936 continuous / 0.972 binary at n =
  800); `reps > 1` adds a small efficiency gain;
  `se_method = "bootstrap"` gives percentile intervals that are no
  better (0.905–0.930) and its `weak_id` flag is not a coverage
  diagnostic; `oe_regular` and the Fieller set handle the near-null
  case, where the ratio is non-regular.
