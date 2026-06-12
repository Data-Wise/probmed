# Build design matrices and free-fit summaries for the parallel MBCO fit.

Generalizes `.mbco_prep` to k parallel mediators. `Dm_list[[j]]` is the
design for mediator j (intercept + its covariates) with the treatment
term removed; `Dy` is the outcome design (intercept + treatment +
outcome covariates) with all mediator terms removed. The b-paths
multiply the columns of `M` (n x k).

## Usage

``` r
.mbco_prep_parallel(extract, x_ref, x_value)
```
