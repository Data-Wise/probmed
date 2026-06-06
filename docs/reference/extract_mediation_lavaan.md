# Extract Mediation Structure from lavaan Objects

Extract mediation structure from SEM models fitted with the `lavaan`
package. Supports models fitted with
[`sem()`](https://rdrr.io/pkg/lavaan/man/sem.html),
[`cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html), or
[`lavaan()`](https://rdrr.io/pkg/lavaan/man/lavaan.html).

## Usage

``` r
extract_mediation_lavaan(
  object,
  treatment,
  mediator,
  outcome = NULL,
  standardized = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted lavaan object.

- treatment:

  Character: name of the treatment variable.

- mediator:

  Character: name of the mediator variable.

- outcome:

  Character: name of the outcome variable. If `NULL`, the function
  attempts to detect it automatically.

- standardized:

  Logical: whether to use standardized estimates. Default is `FALSE`.

- ...:

  Additional arguments (currently unused).

## Value

A `MediationExtract` object.

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("lavaan", quietly = TRUE)) {
    library(lavaan)

    # Define model
    model <- "
    M ~ a*X
    Y ~ b*M + cp*X
  "

    # Simulate data
    X <- rnorm(100)
    M <- 0.5 * X + rnorm(100)
    Y <- 0.4 * M + 0.2 * X + rnorm(100)
    data <- data.frame(X = X, M = M, Y = Y)

    # Fit model
    fit <- sem(model, data = data)

    # Extract
    extract <- extract_mediation(fit, treatment = "X", mediator = "M")

    # Compute P_med
    pmed(extract)
}
} # }
```
