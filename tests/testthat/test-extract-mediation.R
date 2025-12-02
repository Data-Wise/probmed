test_that("extract_mediation works with mediate objects (Linear)", {
    skip_if_not_installed("mediation")

    # Generate data
    data <- generate_mediation_data(n = 100, a = 0.5, b = 0.4, c_prime = 0.2, seed = 123)

    # Fit models
    model_m <- lm(M ~ X, data = data)
    model_y <- lm(Y ~ M + X, data = data)

    # Run mediate
    # We use boot=FALSE for speed, as we only need the object structure
    med_out <- mediation::mediate(model_m, model_y, treat = "X", mediator = "M", boot = FALSE)

    # Extract
    extract <- extract_mediation(med_out)

    # Validate
    expect_true(S7::S7_inherits(extract, MediationExtract))
    expect_equal(extract@treatment, "X")
    expect_equal(extract@mediator, "M")
    expect_equal(extract@outcome, "Y")

    # Check estimates match coefficients
    expect_equal(extract@a_path, coef(model_m)["X"], ignore_attr = TRUE)
    expect_equal(extract@b_path, coef(model_y)["M"], ignore_attr = TRUE)
})

test_that("extract_mediation works with mediate objects (GLM)", {
    skip_if_not_installed("mediation")

    # Generate data with binary outcome
    data <- generate_mediation_data(n = 200, seed = 456)
    # Convert Y to binary
    data$Y_prob <- plogis(0.4 * data$M + 0.2 * data$X)
    data$Y_bin <- rbinom(200, 1, data$Y_prob)

    # Fit models
    model_m <- lm(M ~ X, data = data)
    model_y <- glm(Y_bin ~ M + X, data = data, family = binomial)

    # Run mediate
    med_out <- mediation::mediate(model_m, model_y, treat = "X", mediator = "M", boot = FALSE)

    # Extract
    extract <- extract_mediation(med_out)

    # Validate
    expect_true(S7::S7_inherits(extract, MediationExtract))
    expect_equal(extract@outcome, "Y_bin")

    # Check estimates
    expect_equal(extract@b_path, coef(model_y)["M"], ignore_attr = TRUE)
})

test_that("extract_mediation handles missing models in mediate object", {
    skip_if_not_installed("mediation")

    # Create a fake mediate object without models
    fake_med <- list(treat = "X", mediator = "M")
    class(fake_med) <- "mediate"

    expect_error(
        extract_mediation(fake_med),
        "The mediate object does not contain the underlying models"
    )
})
test_that("extract_mediation handles covariates in mediate objects", {
    skip_if_not_installed("mediation")

    data <- generate_mediation_data(n = 100, n_covariates = 2, seed = 789)

    model_m <- lm(M ~ X + C1 + C2, data = data)
    model_y <- lm(Y ~ M + X + C1 + C2, data = data)

    med_out <- mediation::mediate(model_m, model_y, treat = "X", mediator = "M", boot = FALSE)

    extract <- extract_mediation(med_out)

    expect_true("C1" %in% extract@mediator_predictors)
    expect_true("C2" %in% extract@mediator_predictors)
    expect_true("C1" %in% extract@outcome_predictors)
})

test_that("extract_mediation matches direct lm extraction", {
    skip_if_not_installed("mediation")

    data <- generate_mediation_data(n = 100, seed = 101)

    model_m <- lm(M ~ X, data = data)
    model_y <- lm(Y ~ M + X, data = data)

    med_out <- mediation::mediate(model_m, model_y, treat = "X", mediator = "M", boot = FALSE)

    extract_med <- extract_mediation(med_out)
    extract_lm <- extract_mediation(model_m, model_y = model_y, treatment = "X", mediator = "M")

    expect_equal(extract_med@a_path, extract_lm@a_path)
    expect_equal(extract_med@b_path, extract_lm@b_path)
    expect_equal(extract_med@c_prime, extract_lm@c_prime)
    expect_equal(extract_med@sigma_m, extract_lm@sigma_m)
})

test_that("extract_mediation handles interactions (basic support)", {
    skip_if_not_installed("mediation")

    data <- generate_mediation_data(n = 100, seed = 202)
    # Add interaction effect to Y
    data$Y <- 0.4 * data$M + 0.2 * data$X + 0.1 * data$M * data$X + rnorm(100)

    model_m <- lm(M ~ X, data = data)
    model_y <- lm(Y ~ M * X, data = data) # Interaction M*X

    med_out <- mediation::mediate(model_m, model_y, treat = "X", mediator = "M", boot = FALSE)

    # Currently, probmed's extract_mediation.lm extracts 'M' coefficient as b_path.
    # If there is an interaction, the effect of M depends on X.
    # extract_mediation.lm might extract the main effect of M.
    # This test verifies behavior, even if it's "naive" extraction for now.

    extract <- extract_mediation(med_out)

    # Check if b_path is the main effect of M
    expect_equal(extract@b_path, coef(model_y)["M"], ignore_attr = TRUE)

    # Note: P_med calculation currently assumes no interaction.
    # Future work should handle interactions explicitly.
})
