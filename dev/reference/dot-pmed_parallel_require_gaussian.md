# Require a Gaussian parallel extract (finite residual SDs present).

`ParallelMediationData` carries no `family_*` slots, so the Gaussian
indicator is "all `sigma_mediators` and `sigma_y` present and finite". A
non-Gaussian GLM fit leaves these NA/empty. Caveat: a lavaan-sourced
extract computes a residual SD for any variable under normal theory, so
a binary / ordinal mediator fit that way could slip past this heuristic
– callers must ensure mediators and outcome are genuinely
continuous-Gaussian. The durable fix is family slots on the medfit
ParallelMediationData container.

## Usage

``` r
.pmed_parallel_require_gaussian(extract)
```
