# Model terms that involve a variable other than as its own main effect

Returns the term labels (e.g. `"X:M"`, `"I(M^2)"`, `"log(X)"`) that
mention any of `vars` without being exactly one of them.

## Usage

``` r
.offending_terms(labels, vars)
```
