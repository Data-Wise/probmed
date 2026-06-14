# Deterministic perturbed warm starts for the parallel MBCO optimizers.

Returns `base` followed by `n_extra` jittered copies. The jitter draws
from a fixed local RNG seed whose prior global state is restored on
exit, so the starts (and hence the MBCO interval) stay reproducible and
seed-free for the caller. The constrained surfaces are multimodal /
penalty-walled, so a single warm start can stick below the true
constrained MLE; multistart guards that.

## Usage

``` r
.mbco_parallel_starts(base, n_extra, sd = 0.25)
```
