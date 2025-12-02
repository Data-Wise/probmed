#' Extract Mediation Structure from mediate Objects
#'
#' @description
#' Extract mediation structure from objects produced by the \code{mediation} package.
#' This method extracts the underlying \code{lm} or \code{glm} models and delegates
#' to the corresponding extraction methods.
#'
#' @param object A \code{mediate} object.
#' @param ... Additional arguments passed to the underlying extraction method.
#'
#' @return A \code{MediationExtract} object.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("mediation", quietly = TRUE)) {
#'     library(mediation)
#'
#'     # Generate data
#'     set.seed(123)
#'     n <- 100
#'     data <- data.frame(X = rnorm(n), C = rnorm(n))
#'     data$M <- 0.5 * data$X + 0.3 * data$C + rnorm(n)
#'     data$Y <- 0.4 * data$M + 0.2 * data$X + 0.1 * data$C + rnorm(n)
#'
#'     # Fit models
#'     model_m <- lm(M ~ X + C, data = data)
#'     model_y <- lm(Y ~ M + X + C, data = data)
#'
#'     # Run mediate
#'     med_out <- mediate(model_m, model_y, treat = "X", mediator = "M", boot = FALSE)
#'
#'     # Extract
#'     extract <- extract_mediation(med_out)
#'
#'     # Compute P_med
#'     pmed(extract)
#' }
#' }
#'
#' @export
S7::method(extract_mediation, mediate_class) <- function(object, ...) {
    # Extract underlying models
    model_m <- object$model.m
    model_y <- object$model.y

    if (is.null(model_m) || is.null(model_y)) {
        stop("The mediate object does not contain the underlying models (model.m and model.y).")
    }

    # Extract variable names
    treatment <- object$treat
    mediator <- object$mediator

    # Delegate to existing method
    # We call extract_mediation on model_m, passing model_y and names
    # Note: extract_mediation is a generic, so it will dispatch based on class(model_m)
    extract_mediation(
        model_m,
        model_y = model_y,
        treatment = treatment,
        mediator = mediator,
        ...
    )
}
