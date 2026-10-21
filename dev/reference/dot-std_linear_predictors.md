# Per-row covariate linear predictors at x = 0 and m = 0

`coef_m` / `coef_y` are named by design column (no prefix). Returns the
mediator and outcome linear predictors with the treatment and mediator
contributions removed, one value per analysis row.

## Usage

``` r
.std_linear_predictors(design, coef_m, coef_y, treatment, mediator)
```
