# P_med Result Class

S7 class for P_med computation results.

## Usage

``` r
PmedResult(
  estimate = integer(0),
  ci_lower = integer(0),
  ci_upper = integer(0),
  ci_level = integer(0),
  method = character(0),
  n_boot = NA_integer_,
  boot_estimates = numeric(0),
  ie_estimate = NA_real_,
  ie_ci_lower = NA_real_,
  ie_ci_upper = NA_real_,
  ie_boot_estimates = numeric(0),
  x_ref = integer(0),
  x_value = integer(0),
  source_extract = MediationExtract(),
  converged = logical(0),
  call = NULL
)
```

## Arguments

- estimate:

  Numeric: P_med point estimate

- ci_lower:

  Numeric: lower bound of confidence interval

- ci_upper:

  Numeric: upper bound of confidence interval

- ci_level:

  Numeric: confidence level (e.g., 0.95)

- method:

  Character: inference method used

- n_boot:

  Integer: number of bootstrap samples (NA if not bootstrap)

- boot_estimates:

  Numeric vector: bootstrap distribution (empty if not bootstrap)

- ie_estimate:

  Numeric: Indirect Effect (NIE) point estimate

- ie_ci_lower:

  Numeric: lower bound of NIE confidence interval

- ie_ci_upper:

  Numeric: upper bound of NIE confidence interval

- ie_boot_estimates:

  Numeric vector: bootstrap distribution of NIE

- x_ref:

  Numeric: reference treatment value

- x_value:

  Numeric: treatment value for contrast

- source_extract:

  MediationExtract object: source of the estimates

- converged:

  Logical: did computation converge

- call:

  Call object: original function call
