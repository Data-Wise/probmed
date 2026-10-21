# Convergence-checked nloptr minimizer with one deterministic fallback

Derivative-free BOBYQA first, then Subplex (a robust Nelder-Mead
variant) if BOBYQA does not report success. Returns the shape the
callers use: `par`, `value`, and `convergence` (0 on an nloptr success
status 1-4).

## Usage

``` r
.mbco_optim(start, fn)
```
