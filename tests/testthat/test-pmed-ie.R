test_that("pmed computes Indirect Effect (IE)", {
    # Generate data
    set.seed(123)
    n <- 100
    data <- data.frame(X = rnorm(n))
    data$M <- 0.5 * data$X + rnorm(n)
    data$Y <- 0.4 * data$M + 0.2 * data$X + rnorm(n)

    # Fit models
    model_m <- lm(M ~ X, data = data)
    model_y <- lm(Y ~ M + X, data = data)

    # Extract
    extract <- extract_mediation(model_m, model_y = model_y, treatment = "X", mediator = "M", data = data)

    # 1. Plugin method
    res_plugin <- pmed(extract, method = "plugin")

    # Expected IE = a * b
    a <- coef(model_m)["X"]
    b <- coef(model_y)["M"]
    expected_ie <- as.numeric(a * b)

    expect_equal(res_plugin@ie_estimate, expected_ie)

    # 2. Parametric Bootstrap
    if (requireNamespace("MASS", quietly = TRUE)) {
        res_boot <- pmed(extract, method = "parametric_bootstrap", n_boot = 100, seed = 456)

        # Check IE estimate exists
        expect_true(!is.na(res_boot@ie_estimate))
        expect_true(!is.na(res_boot@ie_ci_lower))
        expect_true(!is.na(res_boot@ie_ci_upper))

        # Check if IE estimate is in reasonable range
        # Bootstrap has randomness - just check sign and rough magnitude
        expect_true(abs(res_boot@ie_estimate - expected_ie) < 1.0)
    }

    # 3. Nonparametric Bootstrap
    res_np <- pmed(extract, method = "nonparametric_bootstrap", n_boot = 50, seed = 789)
    expect_true(!is.na(res_np@ie_estimate))
    expect_true(!is.na(res_np@ie_ci_lower))
    expect_true(!is.na(res_np@ie_ci_upper))
})
