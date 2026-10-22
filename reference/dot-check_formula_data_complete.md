# Reject missing values in the user's data before any model is fitted

Newer medfit extractions keep only the rows a fit used, so a later check
on the extracted data cannot see the missing rows (and the two models
can be fitted to different samples). The check therefore runs on the
data the caller passed, using the variables named in either formula.

## Usage

``` r
.check_formula_data_complete(formula_y, formula_m, data)
```
