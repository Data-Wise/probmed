# Simple P_med Core Computation

Computes the mediation estimand \\P\_{med} = P(Y(x, M(x)) \> Y(x,
M(x^\*))) + \tfrac12 P(Y(x, M(x)) = Y(x, M(x^\*)))\\: the mediator is
drawn under BOTH treatment levels independently while the treatment is
held at `x_value` for both outcomes, so the direct effect `c_prime`
cancels and the result reflects mediation only. The intercepts and
families matter only for non-Gaussian models (they cancel in the
Gaussian linear contrast); their defaults keep Gaussian behaviour
unchanged.

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
  i_m = 0
)
```
