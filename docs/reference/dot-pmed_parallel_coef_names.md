# Names of the structural coefficients (a_j, b_j) in the parallel estimates vector / vcov. Uses the structural aliases `a{j}` / `b{j}`, which BOTH the lm/glm and lavaan medfit extractors expose (with full vcov rows) – unlike the source-specific `m{j}_<tx>` / `y_<mediator_j>` names, which exist only in the lm/glm extractor and would break the parametric bootstrap on a lavaan extract.

Names of the structural coefficients (a_j, b_j) in the parallel
estimates vector / vcov. Uses the structural aliases `a{j}` / `b{j}`,
which BOTH the lm/glm and lavaan medfit extractors expose (with full
vcov rows) – unlike the source-specific `m{j}_<tx>` / `y_<mediator_j>`
names, which exist only in the lm/glm extractor and would break the
parametric bootstrap on a lavaan extract.

## Usage

``` r
.pmed_parallel_coef_names(extract)
```
