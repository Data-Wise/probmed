# Build design matrices and free-fit summaries for the MBCO optimization.

`Dm` is the M-equation design (intercept plus mediator covariates) with
the treatment term `a*X` removed; `Dy` is the Y-equation design
(intercept, treatment, and outcome covariates) with the mediator term
`b*M` removed.

## Usage

``` r
.mbco_prep(extract, x_ref, x_value)
```
