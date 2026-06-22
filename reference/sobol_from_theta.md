# Map corner means to the Sobol variance decomposition

Closed-form map from the four corner means `theta(a, a')` to the Sobol /
functional-ANOVA variance components and the variance-scale proportion
mediated. Used internally by
[`sobol_pmed()`](https://data-wise.github.io/probmed/reference/sobol_pmed.md)
and exported as a helper for computing the analytic (population) truth
from known corner means.

## Usage

``` r
sobol_from_theta(theta, pd = 0.5, pm = 0.5)
```

## Arguments

- theta:

  Numeric length-4: corner means named `"11"`, `"10"`, `"01"`, `"00"`
  (i.e. `theta(a, a')` for `(a, a')` in the four corners).

- pd, pm:

  Numeric: tilt probabilities for the direct and mediated pathways
  (Bernoulli design points); defaults `0.5`.

## Value

A named numeric vector with elements `Vd`, `Vm`, `Vdm`, `VT`,
`Pmed_sobol` (= `Vm / VT`), `S1_med` (= `Vm / VT`), `ST_med` (=
`(Vm + Vdm) / VT`), and `Rint` (the interaction remainder).

## Examples

``` r
th <- c("11" = 1.32, "10" = 0.5, "01" = 0.42, "00" = 0)
sobol_from_theta(th)
#>         Vd         Vm        Vdm         VT Pmed_sobol     S1_med     ST_med 
#>  0.1225000  0.0961000  0.0100000  0.2286000  0.4203850  0.4203850  0.4641295 
#>       Rint 
#>  0.4000000 
```
