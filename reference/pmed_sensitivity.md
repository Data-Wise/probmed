# M–Y Confounding Sensitivity for a Proportion-Mediated Estimand

For a mean-based proportion mediated, this function adds a specified
bias `b` to the indirect effect while holding the total effect fixed. It
reports `(indirect + b) / total` over the bias grid and the bias needed
to make the indirect effect zero or bring the ratio to a user-specified
threshold.

## Usage

``` r
pmed_sensitivity(
  object = NULL,
  indirect = NULL,
  total = NULL,
  bias_grid = NULL,
  threshold = NA_real_,
  evalue = FALSE,
  se = NA_real_,
  scale = NA_character_,
  ...
)
```

## Arguments

- object:

  A
  [GaugePmedResult](https://data-wise.github.io/probmed/reference/GaugePmedResult.md)
  or another object exposing indirect and total effect components.
  Ignored if both `indirect` and `total` are supplied directly.

- indirect, total:

  Numeric: supply the components directly instead of (or to override)
  extraction from `object`.

- bias_grid:

  Numeric: numerator biases `b` to evaluate. Default spans
  `+/- |indirect|` on an 11-point grid.

- threshold:

  Numeric or `NA`: a ratio value whose tipping bias to report (e.g. `0`
  to ask "how much bias makes the share vanish"). Default `NA` (skip).

- evalue:

  Logical: also compute an EValue E-value summary if the package is
  installed (default `FALSE`).

- se:

  Numeric or `NA`: standard error of `indirect`, used only for the
  E-value CI bound when `evalue = TRUE`.

- scale:

  Character: label for the `P_med` scale (free text).

- ...:

  Unused.

## Value

A
[PmedSensitivity](https://data-wise.github.io/probmed/reference/PmedSensitivity.md)
object.

## Details

The additive bias model applies to the mean-based gauge and incremental
shares from
[`ward_residual()`](https://data-wise.github.io/probmed/reference/ward_residual.md)
and
[`incr_pmed()`](https://data-wise.github.io/probmed/reference/incr_pmed.md).
In those frameworks the bias represents a shift in the cross-world mean
when their assumption A4 (no exposure-induced mediator–outcome
confounding) fails.

The tie-adjusted stochastic-superiority probability from
[`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md) and
the manuscript's A2/A5 sensitivity analyses are outside the scope of
this helper. For a Sobol variance share, callers would need to supply a
separately justified perturbation of the variance component.

The optional E-value summary (`evalue = TRUE`) is a coarse RR-oriented
robustness number; it is reported only as a familiar single figure and
is guarded behind the EValue package (Suggests).

## Examples

``` r
## From mean-based effect components:
pmed_sensitivity(indirect = 0.30, total = 0.50, threshold = 0)
#> P_med M-Y confounding sensitivity (A4)
#>   P_med = 0.600  (indirect=0.300 / total=0.500)
#>   Tipping bias (indirect -> 0): b = -0.300
#>   Tipping bias (P_med -> 0.000): b = -0.300
#>   P_med over bias grid [-0.300, 0.300]: [0.000, 1.200]

# \donttest{
## From a Paper 1 fit:
set.seed(1)
n <- 400; C <- rnorm(n)
A <- rbinom(n, 1, plogis(0.5 * C))
M <- 0.5 * A + 0.3 * C + rnorm(n)
Y <- 0.4 * A + 0.6 * M + 0.2 * C + rnorm(n)
fit <- ward_residual(data.frame(A, M, Y, C))
pmed_sensitivity(fit, threshold = 0)
#> P_med M-Y confounding sensitivity (A4)
#>   P_med = 0.514  (indirect=0.338 / total=0.658)
#>   Tipping bias (indirect -> 0): b = -0.338
#>   Tipping bias (P_med -> 0.000): b = -0.338
#>   P_med over bias grid [-0.338, 0.338]: [0.000, 1.028]
# }
```
