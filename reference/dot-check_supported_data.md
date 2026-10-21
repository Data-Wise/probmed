# Data-level checks: missing rows, variable types, binary responses

The nonparametric bootstrap resamples `extract@data`, so rows a model
fit dropped for missingness would come back; and a binomial response
that is not coded 0/1 (proportions, successes/failures) cannot be
simulated as Bernoulli.

## Usage

``` r
.check_supported_data(extract)
```
