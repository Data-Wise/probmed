#' Generate Mediation Data with Known Effects
#'
#' @description
#' Utility function for generating test data with known mediation effects.
#' Useful for validating extraction methods and P_med computations.
#'
#' @param n Sample size (default: 200)
#' @param a Treatment -> Mediator effect (a path, default: 0.5)
#' @param b Mediator -> Outcome effect (b path, default: 0.4)
#' @param c_prime Direct effect of treatment on outcome (c' path, default: 0.2)
#' @param sigma_m Residual standard deviation of mediator model (default: 1)
#' @param sigma_y Residual standard deviation of outcome model (default: 1)
#' @param n_covariates Number of covariates to include (default: 0)
#' @param covariate_effects List with elements 'm' and 'y' specifying covariate
#'   effects on mediator and outcome (default: NULL, uses random effects)
#' @param seed Random seed for reproducibility (default: NULL)
#'
#' @return Data frame with columns X (treatment), M (mediator), Y (outcome),
#'   and C1, C2, ... (covariates if n_covariates > 0). The data frame has
#'   an attribute 'true_effects' containing the true parameter values.
#'
#' @examples
#' # Simple mediation with no covariates
#' data <- generate_mediation_data(n = 100, a = 0.5, b = 0.4, c_prime = 0.2)
#'
#' # With covariates
#' data <- generate_mediation_data(
#'     n = 200,
#'     a = 0.5,
#'     b = 0.4,
#'     c_prime = 0.2,
#'     n_covariates = 2
#' )
#'
#' # Access true effects
#' attr(data, "true_effects")
#'
#' @keywords internal
generate_mediation_data <- function(n = 200,
                                    a = 0.5,
                                    b = 0.4,
                                    c_prime = 0.2,
                                    sigma_m = 1,
                                    sigma_y = 1,
                                    n_covariates = 0,
                                    covariate_effects = NULL,
                                    seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # Treatment (standardized)
    X <- stats::rnorm(n)

    # Covariates
    C_mat <- NULL
    gamma_m <- NULL # Covariate effects on M
    gamma_y <- NULL # Covariate effects on Y

    if (n_covariates > 0) {
        C_mat <- matrix(stats::rnorm(n * n_covariates), ncol = n_covariates)
        colnames(C_mat) <- paste0("C", seq_len(n_covariates))

        # Covariate effects
        if (is.null(covariate_effects)) {
            # Random small effects
            gamma_m <- stats::rnorm(n_covariates, mean = 0, sd = 0.3)
            gamma_y <- stats::rnorm(n_covariates, mean = 0, sd = 0.3)
        } else {
            gamma_m <- covariate_effects$m
            gamma_y <- covariate_effects$y
        }
    }

    # Mediator: M = a*X + gamma_m'*C + epsilon_m
    M <- a * X
    if (n_covariates > 0) {
        M <- M + as.vector(C_mat %*% gamma_m)
    }
    M <- M + stats::rnorm(n, sd = sigma_m)

    # Outcome: Y = b*M + c'*X + gamma_y'*C + epsilon_y
    Y <- b * M + c_prime * X
    if (n_covariates > 0) {
        Y <- Y + as.vector(C_mat %*% gamma_y)
    }
    Y <- Y + stats::rnorm(n, sd = sigma_y)

    # Construct data frame
    data <- data.frame(X = X, M = M, Y = Y)
    if (n_covariates > 0) {
        data <- cbind(data, as.data.frame(C_mat))
    }

    # Store true effects as attribute
    true_effects <- list(
        a = a,
        b = b,
        c_prime = c_prime,
        sigma_m = sigma_m,
        sigma_y = sigma_y,
        indirect_effect = a * b,
        total_effect = a * b + c_prime
    )

    if (n_covariates > 0) {
        true_effects$gamma_m <- gamma_m
        true_effects$gamma_y <- gamma_y
    }

    attr(data, "true_effects") <- true_effects

    return(data)
}


#' Generate Binary Mediation Data
#'
#' @description
#' Generate data with binary mediator and/or outcome for testing GLM extraction.
#'
#' @param n Sample size
#' @param a Treatment effect on mediator (on logit scale if binary)
#' @param b Mediator effect on outcome (on logit scale if binary)
#' @param c_prime Direct effect (on logit scale if binary outcome)
#' @param binary_mediator Logical: should mediator be binary?
#' @param binary_outcome Logical: should outcome be binary?
#' @param seed Random seed
#'
#' @return Data frame with binary variables as appropriate
#'
#' @keywords internal
generate_binary_mediation_data <- function(n = 200,
                                           a = 0.5,
                                           b = 0.5,
                                           c_prime = 0.2,
                                           binary_mediator = TRUE,
                                           binary_outcome = TRUE,
                                           seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # Treatment
    X <- stats::rnorm(n)

    # Mediator
    if (binary_mediator) {
        # Logit model: logit(P(M=1)) = a*X
        p_m <- plogis(a * X)
        M <- stats::rbinom(n, 1, p_m)
    } else {
        # Continuous
        M <- a * X + stats::rnorm(n)
    }

    # Outcome
    if (binary_outcome) {
        # Logit model: logit(P(Y=1)) = b*M + c'*X
        p_y <- plogis(b * M + c_prime * X)
        Y <- stats::rbinom(n, 1, p_y)
    } else {
        # Continuous
        Y <- b * M + c_prime * X + stats::rnorm(n)
    }

    data <- data.frame(X = X, M = M, Y = Y)

    attr(data, "true_effects") <- list(
        a = a,
        b = b,
        c_prime = c_prime,
        binary_mediator = binary_mediator,
        binary_outcome = binary_outcome
    )

    return(data)
}
