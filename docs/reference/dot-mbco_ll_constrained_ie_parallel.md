# Maximized log-likelihood under the constraint sum(a\*b) = ie0 (parallel).

Pins the last b-path via `b_k = (ie0 - sum_{j<k} a_j b_j) / a_k`; Vy is
free (MLE = SSE_Y / n). Optimizes over
`(a_1..k, log Vm_1..k, b_1..k-1)`. Because the b_k = .../a_k
substitution makes the profile multimodal away from the point estimate,
warm-starts from two basins (as the single-mediator IE fit).

## Usage

``` r
.mbco_ll_constrained_ie_parallel(prep, ie0)
```
