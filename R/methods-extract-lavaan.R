#' Extract Mediation Structure from lavaan Objects
#'
#' @description
#' Extract mediation structure from SEM models fitted with the \code{lavaan} package.
#' Supports models fitted with \code{sem()}, \code{cfa()}, or \code{lavaan()}.
#'
#' @param object A fitted lavaan object.
#' @param treatment Character: name of the treatment variable.
#' @param mediator Character: name of the mediator variable.
#' @param outcome Character: name of the outcome variable. If \code{NULL},
#'   the function attempts to detect it automatically.
#' @param standardized Logical: whether to use standardized estimates.
#'   Default is \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{MediationExtract} object.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#'     library(lavaan)
#'
#'     # Define model
#'     model <- "
#'     M ~ a*X
#'     Y ~ b*M + cp*X
#'   "
#'
#'     # Simulate data
#'     X <- rnorm(100)
#'     M <- 0.5 * X + rnorm(100)
#'     Y <- 0.4 * M + 0.2 * X + rnorm(100)
#'     data <- data.frame(X = X, M = M, Y = Y)
#'
#'     # Fit model
#'     fit <- sem(model, data = data)
#'
#'     # Extract
#'     extract <- extract_mediation(fit, treatment = "X", mediator = "M")
#'
#'     # Compute P_med
#'     pmed(extract)
#' }
#' }
#' @keywords internal
extract_mediation_lavaan <- function(object,
                                     treatment,
                                     mediator,
                                     outcome = NULL,
                                     standardized = FALSE,
                                     ...) {
    # 1. Check convergence
    if (!lavaan::lavInspect(object, "converged")) {
        warning("lavaan model did not converge. Results may be unreliable.")
    }

    # 2. Get parameter estimates
    params <- lavaan::parameterEstimates(object, standardized = standardized)

    # Determine which column to use for estimates
    est_col <- if (standardized) "std.all" else "est"

    # 3. Validate variables and auto-detect outcome

    # Check if treatment and mediator exist in the model
    all_vars <- lavaan::lavNames(object)
    if (!treatment %in% all_vars) {
        stop("Treatment variable '", treatment, "' not found in the model.")
    }
    if (!mediator %in% all_vars) {
        stop("Mediator variable '", mediator, "' not found in the model.")
    }

    # Auto-detect outcome if needed
    if (is.null(outcome)) {
        # Outcome should be endogenous (lhs of ~) and have both M and X as predictors
        # Find variables where M is a predictor
        has_mediator <- params[params$op == "~" & params$rhs == mediator, "lhs"]
        # Find variables where X is a predictor
        has_treatment <- params[params$op == "~" & params$rhs == treatment, "lhs"]

        potential_outcomes <- intersect(has_mediator, has_treatment)
        # Exclude mediator itself (if M ~ X is in model)
        potential_outcomes <- setdiff(potential_outcomes, mediator)

        if (length(potential_outcomes) == 0) {
            stop("Could not auto-detect outcome variable. Please specify 'outcome'.")
        } else if (length(potential_outcomes) > 1) {
            stop(
                "Multiple potential outcomes detected (",
                paste(potential_outcomes, collapse = ", "),
                "). Please specify 'outcome'."
            )
        }
        outcome <- potential_outcomes
    } else {
        if (!outcome %in% all_vars) {
            stop("Outcome variable '", outcome, "' not found in the model.")
        }
    }

    # 4. Extract path coefficients

    # a path: M ~ X
    a_row <- params[params$lhs == mediator & params$rhs == treatment & params$op == "~", ]
    if (nrow(a_row) == 0) {
        stop("Path ", mediator, " ~ ", treatment, " (a path) not found in model.")
    }
    a_path <- a_row[[est_col]]

    # b path: Y ~ M
    b_row <- params[params$lhs == outcome & params$rhs == mediator & params$op == "~", ]
    if (nrow(b_row) == 0) {
        stop("Path ", outcome, " ~ ", mediator, " (b path) not found in model.")
    }
    b_path <- b_row[[est_col]]

    # c' path: Y ~ X
    cp_row <- params[params$lhs == outcome & params$rhs == treatment & params$op == "~", ]
    if (nrow(cp_row) == 0) {
        stop("Path ", outcome, " ~ ", treatment, " (c' path) not found in model.")
    }
    c_prime <- cp_row[[est_col]]

    # 5. Extract Variance-Covariance Matrix

    # We need the vcov for the coefficients a, b, and c'
    # lavaan::vcov() returns the full matrix. We need to subset it.
    # We use the parameter table to find the indices of the free parameters.

    partable <- lavaan::parTable(object)

    # Helper to find free parameter index
    get_param_idx <- function(lhs, op, rhs) {
        idx <- which(partable$lhs == lhs & partable$op == op & partable$rhs == rhs)
        if (length(idx) == 0) {
            return(NA)
        }
        # Check if it is a free parameter
        if (partable$free[idx] == 0) {
            return(NA)
        } # Fixed parameter
        return(partable$free[idx])
    }

    a_idx <- get_param_idx(mediator, "~", treatment)
    b_idx <- get_param_idx(outcome, "~", mediator)
    cp_idx <- get_param_idx(outcome, "~", treatment)

    # If any path is fixed, its variance is 0 (and covs are 0)
    # We construct the 3x3 vcov matrix manually
    vcov_full <- lavaan::vcov(object)
    vcov_subset <- matrix(0, 3, 3)
    rownames(vcov_subset) <- colnames(vcov_subset) <- c("a", "b", "cp")

    indices <- c(a = a_idx, b = b_idx, cp = cp_idx)

    # Fill the matrix
    for (i in 1:3) {
        for (j in 1:3) {
            idx_i <- indices[i]
            idx_j <- indices[j]

            if (!is.na(idx_i) && !is.na(idx_j)) {
                vcov_subset[i, j] <- vcov_full[idx_i, idx_j]
            }
        }
    }

    # Note: If standardized=TRUE, the vcov of standardized estimates is different.
    # lavaan::vcov() returns vcov of unstandardized estimates.
    # Getting vcov of standardized estimates is complex (delta method).
    # For now, we warn if standardized=TRUE and vcov is requested for bootstrap.
    if (standardized) {
        warning(
            "Variance-covariance matrix corresponds to unstandardized estimates. ",
            "Standard errors for standardized estimates may be incorrect if used directly."
        )
        # Ideally, we should use lavaan::standardizedSolution(..., se = TRUE) but that gives SEs, not full vcov easily without more work.
        # For P_med plugin, vcov is not used. For bootstrap, we might need to rethink.
        # However, the user might just want the point estimate.
    }

    # 6. Extract Residual Variances (Sigma)

    # Sigma M: M ~~ M
    sigma_m_row <- params[params$lhs == mediator & params$rhs == mediator & params$op == "~~", ]
    sigma_m <- if (nrow(sigma_m_row) > 0) sqrt(sigma_m_row[[est_col]]) else NA_real_

    # Sigma Y: Y ~~ Y
    sigma_y_row <- params[params$lhs == outcome & params$rhs == outcome & params$op == "~~", ]
    sigma_y <- if (nrow(sigma_y_row) > 0) sqrt(sigma_y_row[[est_col]]) else NA_real_

    # 7. Extract Data
    data <- tryCatch(
        {
            d <- lavaan::lavInspect(object, "data")
            if (is.list(d)) {
                # Multi-group: take first group and warn
                warning("Multi-group models not fully supported. Using data from the first group.")
                as.data.frame(d[[1]])
            } else {
                as.data.frame(d)
            }
        },
        error = function(e) {
            warning("Could not extract data from lavaan object: ", e$message)
            NULL
        }
    )

    # 8. Predictors
    mediator_predictors <- params[params$lhs == mediator & params$op == "~", "rhs"]
    outcome_predictors <- params[params$lhs == outcome & params$op == "~", "rhs"]

    # 9. Construct Object
    MediationExtract(
        estimates = c(a = a_path, b = b_path, cp = c_prime),
        mediator_predictors = unique(mediator_predictors),
        outcome_predictors = unique(outcome_predictors),
        a_path = a_path,
        b_path = b_path,
        c_prime = c_prime,
        vcov = vcov_subset,
        sigma_m = sigma_m,
        sigma_y = sigma_y,
        data = data,
        n_obs = lavaan::lavInspect(object, "nobs"),
        source_package = "lavaan",
        converged = lavaan::lavInspect(object, "converged"),
        treatment = treatment,
        mediator = mediator,
        outcome = outcome
    )
}
