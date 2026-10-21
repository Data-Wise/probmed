# Deterministic separation check for a binary-response design

Complete or quasi-complete separation exists exactly when some nonzero
`beta` satisfies `s_i * x_i' beta >= 0` for every row, with
`s_i = 2 y_i - 1`, and strict inequality for at least one row (Albert
and Anderson, 1984). This solves the linear program

## Usage

``` r
.detect_separation(X, y, tol = 1e-09)
```

## Arguments

- X:

  Numeric design matrix (including any intercept column).

- y:

  0/1 response.

- tol:

  Optimum above `tol * max(1, sum(abs(X)))` counts as separation.

## Value

`"separated"`, `"ok"`, or `"unchecked"` when `lpSolve` is not installed
or the LP fails.

## Details

maximize sum_i s_i x_i' beta subject to s_i x_i' beta \>= 0 for all i,
-1 \<= beta_j \<= 1,

(Konis, 2007) with `lpSolve`, splitting `beta = u - v`, `u, v >= 0`. A
positive optimum means separation; zero means none. Fitted probabilities
near 0/1 or glm warnings are not used: they are not proof either way.
