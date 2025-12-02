# Extract Mediation Structure from mediate Objects

Extract mediation structure from objects produced by the `mediation`
package. This method extracts the underlying `lm` or `glm` models and
delegates to the corresponding extraction methods.

## Arguments

- object:

  A `mediate` object.

- ...:

  Additional arguments passed to the underlying extraction method.

## Value

A `MediationExtract` object.

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("mediation", quietly = TRUE)) {
    library(mediation)

    # Generate data
    set.seed(123)
    n <- 100
    data <- data.frame(X = rnorm(n), C = rnorm(n))
    data$M <- 0.5 * data$X + 0.3 * data$C + rnorm(n)
    data$Y <- 0.4 * data$M + 0.2 * data$X + 0.1 * data$C + rnorm(n)

    # Fit models
    model_m <- lm(M ~ X + C, data = data)
    model_y <- lm(Y ~ M + X + C, data = data)

    # Run mediate
    med_out <- mediate(model_m, model_y, treat = "X", mediator = "M", boot = FALSE)

    # Extract
    extract <- extract_mediation(med_out)

    # Compute P_med
    pmed(extract)
}
} # }
```
