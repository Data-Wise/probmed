# Profile LR statistic for H0: P_med = p (single-mediator Gaussian MBCO)

`-2 * (ll_constrained - ll_free)`; the quantity the interval inverts.
`NA` when the constrained profile fails. Used by the composite-null
calibration harness (`inst/sim/calibration/calib_mbco_null.R`).

## Usage

``` r
.mbco_lr_at(prep, p)
```
