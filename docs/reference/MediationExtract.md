# Mediation Extract Base Class

S7 class for extracted mediation structures from fitted models.

## Usage

``` r
MediationExtract(
  estimates = integer(0),
  mediator_predictors = character(0),
  outcome_predictors = character(0),
  a_path = integer(0),
  b_path = integer(0),
  c_prime = integer(0),
  vcov = NULL,
  sigma_m = integer(0),
  sigma_y = integer(0),
  data = NULL,
  n_obs = integer(0),
  source_package = character(0),
  converged = logical(0),
  treatment = character(0),
  mediator = character(0),
  outcome = character(0)
)
```

## Arguments

- estimates:

  Numeric vector of all model parameter estimates

- mediator_predictors:

  Character vector of mediator model predictor names

- outcome_predictors:

  Character vector of outcome model predictor names

- a_path:

  Numeric: treatment effect on mediator (a path)

- b_path:

  Numeric: mediator effect on outcome (b path)

- c_prime:

  Numeric: direct effect of treatment on outcome (c' path)

- vcov:

  Variance-covariance matrix of parameter estimates

- sigma_m:

  Numeric: residual standard deviation of mediator model

- sigma_y:

  Numeric: residual standard deviation of outcome model

- data:

  Data frame used for model fitting

- n_obs:

  Integer: number of observations

- source_package:

  Character: name of package that created the models

- converged:

  Logical: did the model(s) converge

- treatment:

  Character: name of treatment variable

- mediator:

  Character: name of mediator variable

- outcome:

  Character: name of outcome variable
