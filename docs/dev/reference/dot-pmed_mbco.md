# MBCO Confidence Interval for P_med (Gaussian single-mediator model)

Implements the Model-Based Constrained Optimization interval (Tofighi &
Kelley, 2020) by likelihood-ratio inversion. Gaussian outcome AND
mediator only; covariates in either equation are supported as profiled
nuisance parameters. Deterministic (no resampling).

## Usage

``` r
.pmed_mbco(extract, x_ref, x_value, ci_level = 0.95, ...)
```
