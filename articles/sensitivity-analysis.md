# Sensitivity analysis for unmeasured M-Y confounding

Randomizing the treatment does not remove mediator-outcome confounding.
[`pmed_sensitivity()`](https://data-wise.github.io/probmed/reference/pmed_sensitivity.md)
tracks one assumption, **A4** (no exposure-induced mediator-outcome
confounding). It biases the cross-world mean in the *indirect* component
and leaves the *total* effect, which is identified from observable
means, unaffected. The function adds a stated numerator bias $`b`$ and
reports the resulting share $`(\text{indirect} + b) / \text{total}`$
over a grid.

**Scope.** The additive-offset model is exact for the mean-based shares:
the gauge
([`ward_residual()`](https://data-wise.github.io/probmed/reference/ward_residual.md))
and incremental
([`incr_pmed()`](https://data-wise.github.io/probmed/reference/incr_pmed.md))
estimands. It is not the right object for the Sobol variance share, and
it is not applied to the core
[`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md)
probability. Do not read its grid as a corrected estimate.

## From components

If you have an indirect and a total effect from any source:

``` r

library(probmed)
pmed_sensitivity(indirect = 0.30, total = 0.50, threshold = 0, evalue = TRUE)
```

The printed result gives the tipping bias that drives the indirect
effect to zero and the bias that drives the share to `threshold`. With
`evalue = TRUE` it adds a coarse E-value, which needs the suggested
package EValue; `se` (the standard error of `indirect`) is used only for
the E-value’s interval bound.

## From a fitted object

``` r

fit <- ward_residual(data.frame(A, M, Y, C))
pmed_sensitivity(fit, threshold = 0)
```

The default grid spans plus and minus the absolute indirect effect in 11
points; pass `bias_grid` to change it. The share can leave \[0, 1\] on
the grid. That signals an implausibly large bias, not a valid share.

## Incremental version

[`incr_sensitivity()`](https://data-wise.github.io/probmed/reference/incr_sensitivity.md)
gives the same tipping analysis per shift value delta for
[`incr_pmed()`](https://data-wise.github.io/probmed/reference/incr_pmed.md).
See [Incremental
P_med](https://data-wise.github.io/probmed/articles/incremental-pmed.md).

## Reading the result

1.  Pick a threshold that matters substantively (zero is the usual one).
2.  Compare the tipping bias to the size of covariate effects you can
    name. A tipping bias larger than any measured confounder is evidence
    of robustness; a smaller one is not.
3.  Report the E-value, if used, as a coarse summary, not as a test.

This analysis addresses M-Y confounding only. It does not address
treatment-mediator confounding, measurement error, or the boundary
behavior described in [Modern
estimands](https://data-wise.github.io/probmed/articles/modern-estimands.html#boundary-behavior).
