# Invert an LR test: find where excess(t) = LR(t) - crit crosses zero.

Steps outward from `est` in increments of `step` until the statistic
exceeds the cutoff, then locates the crossing with uniroot
(resolution-independent). `domain = c(lo, hi)` bounds the search; an
endpoint pinned at a bound means the interval is open there. The
step-out loop is sized to span the whole domain (`max_k` steps), so a
crossing anywhere in `[lo, hi]` is reached – a fixed iteration cap could
otherwise stop short of a wide domain and return a non-root interior
point.

## Usage

``` r
.mbco_invert(est, excess, domain, step)
```

## Value

c(lower = ., upper = .)
