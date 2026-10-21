# Simple P_med Core Computation

Computes the mediation estimand \\P\_{med} = P(Y(x, M(x)) \> Y(x,
M(x^\*))) + \tfrac12 P(Y(x, M(x)) = Y(x, M(x^\*)))\\: the mediator is
drawn under BOTH treatment levels independently while the treatment is
held at `x_value` for both outcomes, so the direct effect `c_prime`
cancels and the result reflects mediation only. The intercepts and
families matter only for non-Gaussian models (they cancel in the
Gaussian linear contrast); their defaults keep Gaussian behaviour
unchanged. `lp_m0` / `lp_y0`, when given, are per-row covariate linear
predictors (see
[`.std_linear_predictors()`](https://data-wise.github.io/probmed/reference/dot-std_linear_predictors.md));
each pair then draws one analysis row and uses it in both arms,
standardizing over the observed covariates. `u`, when given, is an
`n_sim` x 5 matrix of uniforms (see
[`.crn_uniforms()`](https://data-wise.github.io/probmed/reference/dot-crn_uniforms.md))
used in place of fresh RNG draws via the inverse CDF. Reusing one `u`
across bootstrap draws (common random numbers) makes their Monte Carlo
error move together, so it no longer inflates the bootstrap spread.

## Usage

``` r
.pmed_core_simple(
  a,
  b,
  c_prime,
  x_ref,
  x_value,
  sigma_y = 1,
  sigma_m = 1,
  n_sim = 10000,
  family_y = NULL,
  family_m = NULL,
  i_y = 0,
  i_m = 0,
  lp_m0 = NULL,
  lp_y0 = NULL,
  u = NULL
)
```
