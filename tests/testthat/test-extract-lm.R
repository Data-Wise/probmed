test_that("extract_mediation validates required arguments", {
    # Generate test data
    data <- generate_mediation_data(n = 100, seed = 123)
    fit_m <- lm(M ~ X, data = data)
    fit_y <- lm(Y ~ X + M, data = data)

    # Missing model_y
    expect_error(
        extract_mediation(fit_m, treatment = "X", mediator = "M"),
        "provide both models"
    )
})

test_that("extract_mediation validates treatment variable", {
    data <- generate_mediation_data(n = 100, seed = 123)
    fit_m <- lm(M ~ X, data = data)
    fit_y <- lm(Y ~ X + M, data = data)

    # Treatment not in mediator model
    expect_error(
        extract_mediation(fit_m,
            model_y = fit_y,
            treatment = "NonExistent", mediator = "M"
        ),
        "Treatment variable 'NonExistent' not found in mediator model"
    )

    # Treatment not in outcome model
    fit_y_bad <- lm(Y ~ M, data = data) # Missing X
    expect_error(
        extract_mediation(fit_m,
            model_y = fit_y_bad,
            treatment = "X", mediator = "M"
        ),
        "Treatment variable 'X' not found in outcome model"
    )
})

test_that("extract_mediation validates mediator variable", {
    data <- generate_mediation_data(n = 100, seed = 123)
    fit_m <- lm(M ~ X, data = data)
    fit_y <- lm(Y ~ X + M, data = data)

    # Mediator not in outcome model
    expect_error(
        extract_mediation(fit_m,
            model_y = fit_y,
            treatment = "X", mediator = "NonExistent"
        ),
        "Mediator variable 'NonExistent' not found in outcome model"
    )
})

test_that("extract_mediation works with covariates", {
    data <- generate_mediation_data(n = 150, n_covariates = 2, seed = 456)
    fit_m <- lm(M ~ X + C1 + C2, data = data)
    fit_y <- lm(Y ~ X + M + C1 + C2, data = data)

    extract <- extract_mediation(
        fit_m,
        model_y = fit_y,
        treatment = "X",
        mediator = "M",
        data = data
    )

    expect_true(S7::S7_inherits(extract, MediationExtract))
    expect_equal(extract@treatment, "X")
    expect_equal(extract@mediator, "M")
    expect_length(extract@mediator_predictors, 4) # Intercept + X + C1 + C2
    expect_length(extract@outcome_predictors, 5) # Intercept + X + M + C1 + C2
})

test_that("extract_mediation handles binary mediator (GLM)", {
    data <- generate_binary_mediation_data(
        n = 200,
        binary_mediator = TRUE,
        binary_outcome = FALSE,
        seed = 789
    )

    fit_m <- glm(M ~ X, family = binomial(), data = data)
    fit_y <- lm(Y ~ X + M, data = data)

    extract <- extract_mediation(
        fit_m,
        model_y = fit_y,
        treatment = "X",
        mediator = "M",
        data = data
    )

    # Binary mediator should have NA sigma_m
    expect_true(is.na(extract@sigma_m))
    # Gaussian outcome should have valid sigma_y
    expect_false(is.na(extract@sigma_y))
    expect_gt(extract@sigma_y, 0)
})

test_that("extract_mediation handles binary outcome (GLM)", {
    data <- generate_binary_mediation_data(
        n = 200,
        binary_mediator = FALSE,
        binary_outcome = TRUE,
        seed = 101
    )

    fit_m <- lm(M ~ X, data = data)
    fit_y <- glm(Y ~ X + M, family = binomial(), data = data)

    extract <- extract_mediation(
        fit_m,
        model_y = fit_y,
        treatment = "X",
        mediator = "M",
        data = data
    )

    # Gaussian mediator should have valid sigma_m
    expect_false(is.na(extract@sigma_m))
    expect_gt(extract@sigma_m, 0)
    # Binary outcome should have NA sigma_y
    expect_true(is.na(extract@sigma_y))
})

test_that("extract_mediation handles both binary mediator and outcome", {
    data <- generate_binary_mediation_data(
        n = 200,
        binary_mediator = TRUE,
        binary_outcome = TRUE,
        seed = 202
    )

    fit_m <- glm(M ~ X, family = binomial(), data = data)
    fit_y <- glm(Y ~ X + M, family = binomial(), data = data)

    extract <- extract_mediation(
        fit_m,
        model_y = fit_y,
        treatment = "X",
        mediator = "M",
        data = data
    )

    # Both should have NA sigmas
    expect_true(is.na(extract@sigma_m))
    expect_true(is.na(extract@sigma_y))
    expect_equal(extract@source_package, "stats")
})

test_that("extract_mediation warns when data unavailable", {
    data <- generate_mediation_data(n = 100, seed = 303)

    # Fit models without storing model frame
    fit_m <- lm(M ~ X, data = data, model = FALSE)
    fit_y <- lm(Y ~ X + M, data = data, model = FALSE)

    expect_warning(
        extract <- extract_mediation(
            fit_m,
            model_y = fit_y,
            treatment = "X",
            mediator = "M"
        ),
        "No data available in model object"
    )

    # Extract should still work, but data should be NULL
    expect_null(extract@data)
})

test_that("extract_mediation uses provided data when available", {
    data <- generate_mediation_data(n = 100, seed = 404)

    fit_m <- lm(M ~ X, data = data, model = FALSE)
    fit_y <- lm(Y ~ X + M, data = data, model = FALSE)

    # Provide data explicitly
    extract <- extract_mediation(
        fit_m,
        model_y = fit_y,
        treatment = "X",
        mediator = "M",
        data = data
    )

    expect_false(is.null(extract@data))
    expect_equal(nrow(extract@data), 100)
})

test_that("extract_mediation produces correct vcov structure", {
    data <- generate_mediation_data(n = 150, seed = 505)
    fit_m <- lm(M ~ X, data = data)
    fit_y <- lm(Y ~ X + M, data = data)

    extract <- extract_mediation(
        fit_m,
        model_y = fit_y,
        treatment = "X",
        mediator = "M",
        data = data
    )

    # vcov should be block diagonal
    n_m <- length(extract@mediator_predictors)
    n_y <- length(extract@outcome_predictors)

    expect_equal(dim(extract@vcov), c(n_m + n_y, n_m + n_y))

    # Check block diagonal structure (off-diagonal blocks should be zero)
    upper_right <- extract@vcov[1:n_m, (n_m + 1):(n_m + n_y)]
    lower_left <- extract@vcov[(n_m + 1):(n_m + n_y), 1:n_m]

    expect_true(all(upper_right == 0))
    expect_true(all(lower_left == 0))
})

test_that("extract_mediation extracts correct path coefficients", {
    # Use known effects
    data <- generate_mediation_data(
        n = 500, # Large sample for accuracy
        a = 0.6,
        b = 0.5,
        c_prime = 0.3,
        seed = 606
    )

    fit_m <- lm(M ~ X, data = data)
    fit_y <- lm(Y ~ X + M, data = data)

    extract <- extract_mediation(
        fit_m,
        model_y = fit_y,
        treatment = "X",
        mediator = "M",
        data = data
    )

    # Check that extracted paths are close to true values
    # (with some tolerance for sampling variability)
    expect_equal(extract@a_path, 0.6, tolerance = 0.2)
    expect_equal(extract@b_path, 0.5, tolerance = 0.2)
    expect_equal(extract@c_prime, 0.3, tolerance = 0.2)
})

test_that("extract_mediation handles non-convergence gracefully", {
    # Create problematic data (perfect collinearity)
    data <- data.frame(
        X = 1:50,
        M = 1:50, # Perfect correlation with X
        Y = rnorm(50)
    )

    # This will produce warnings about singularity - suppress them
    suppressWarnings({
        fit_m <- lm(M ~ X, data = data)
        fit_y <- lm(Y ~ X + M, data = data)

        # Should still extract (lm doesn't have converged flag, but glm does)
        extract <- extract_mediation(
            fit_m,
            model_y = fit_y,
            treatment = "X",
            mediator = "M",
            data = data
        )
    })

    expect_true(S7::S7_inherits(extract, MediationExtract))
})

test_that("MediationExtract validator catches invalid inputs", {
    # a_path must be scalar
    expect_error(
        MediationExtract(
            estimates = c(1, 2, 3),
            mediator_predictors = c("(Intercept)", "X"),
            outcome_predictors = c("(Intercept)", "X", "M"),
            a_path = c(0.5, 0.6), # Invalid: vector
            b_path = 0.4,
            c_prime = 0.2,
            vcov = matrix(0, 5, 5),
            sigma_m = 1,
            sigma_y = 1,
            data = NULL,
            n_obs = 100L,
            source_package = "stats",
            converged = TRUE,
            treatment = "X",
            mediator = "M",
            outcome = "Y"
        ),
        "a_path must be a single numeric value"
    )

    # n_obs must be positive
    expect_error(
        MediationExtract(
            estimates = c(1, 2, 3),
            mediator_predictors = c("(Intercept)", "X"),
            outcome_predictors = c("(Intercept)", "X", "M"),
            a_path = 0.5,
            b_path = 0.4,
            c_prime = 0.2,
            vcov = matrix(0, 5, 5),
            sigma_m = 1,
            sigma_y = 1,
            data = NULL,
            n_obs = 0L, # Invalid
            source_package = "stats",
            converged = TRUE,
            treatment = "X",
            mediator = "M",
            outcome = "Y"
        ),
        "n_obs must be positive"
    )
})
