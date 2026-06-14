# Closed-form joint P_med for parallel mediators (Gaussian).

For k parallel mediators with treatment fixed at the treated level for
both potential outcomes,
`Y_t - Y_c ~ N(delta * sum(a*b), 2*sum(b^2*Vm) + 2*Vy)`, so the joint
P_med is `Phi(delta*sum(a*b) / sqrt(2*sum(b^2*Vm) + 2*Vy))`. Recovers
the single-mediator formula at k = 1.

## Usage

``` r
.pmed_parallel_closed(a_vec, b_vec, sigma_m_vec, sigma_y, delta)
```

## Arguments

- a_vec, b_vec:

  Length-k path vectors (X -\> Mj, Mj -\> Y).

- sigma_m_vec:

  Length-k mediator residual SDs.

- sigma_y:

  Outcome residual SD (scalar).

- delta:

  Treatment contrast `x_value - x_ref`.

## Value

Scalar joint P_med in `[0, 1]`.
