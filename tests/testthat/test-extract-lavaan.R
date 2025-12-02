test_that("extract_mediation works with simple lavaan model", {
  skip_if_not_installed("lavaan")

  # Generate data
  data <- generate_mediation_data(n = 200, a = 0.5, b = 0.4, c_prime = 0.2, seed = 123)

  # Fit lavaan model
  model <- "
    M ~ a*X
    Y ~ b*M + cp*X
  "
  fit <- lavaan::sem(model, data = data)

  # Extract
  extract <- extract_mediation(fit, treatment = "X", mediator = "M")

  # Validate
  expect_true(S7::S7_inherits(extract, MediationExtract))
  expect_equal(extract@source_package, "lavaan")
  expect_equal(extract@treatment, "X")
  expect_equal(extract@mediator, "M")
  expect_equal(extract@outcome, "Y") # Auto-detected

  # Check estimates
  expect_equal(extract@a_path, 0.5, tolerance = 0.2)
  expect_equal(extract@b_path, 0.4, tolerance = 0.2)
  expect_equal(extract@c_prime, 0.2, tolerance = 0.2)

  # Check sigma
  expect_true(!is.na(extract@sigma_m))
  expect_true(!is.na(extract@sigma_y))
  expect_equal(extract@sigma_m, 1, tolerance = 0.2)
  expect_equal(extract@sigma_y, 1, tolerance = 0.2)
})

test_that("extract_mediation handles lavaan models with covariates", {
  skip_if_not_installed("lavaan")

  data <- generate_mediation_data(n = 200, n_covariates = 2, seed = 456)

  model <- "
    M ~ a*X + C1 + C2
    Y ~ b*M + cp*X + C1 + C2
  "
  fit <- lavaan::sem(model, data = data)

  extract <- extract_mediation(fit, treatment = "X", mediator = "M")

  expect_true("C1" %in% extract@mediator_predictors)
  expect_true("C2" %in% extract@mediator_predictors)
  expect_true("C1" %in% extract@outcome_predictors)
  expect_true("C2" %in% extract@outcome_predictors)
})

test_that("extract_mediation works without parameter labels", {
  skip_if_not_installed("lavaan")

  data <- generate_mediation_data(n = 200, seed = 789)

  # No labels
  model <- "
    M ~ X
    Y ~ M + X
  "
  fit <- lavaan::sem(model, data = data)

  extract <- extract_mediation(fit, treatment = "X", mediator = "M")

  expect_true(S7::S7_inherits(extract, MediationExtract))
  expect_equal(extract@outcome, "Y")
})

test_that("extract_mediation handles standardized estimates", {
  skip_if_not_installed("lavaan")

  data <- generate_mediation_data(n = 200, seed = 101)

  model <- "
    M ~ X
    Y ~ M + X
  "
  fit <- lavaan::sem(model, data = data)

  # Warning expected about vcov
  expect_warning(
    extract_std <- extract_mediation(
      fit,
      treatment = "X",
      mediator = "M",
      standardized = TRUE
    ),
    "Variance-covariance matrix corresponds to unstandardized estimates"
  )

  # Standardized estimates should be between -1 and 1
  expect_gte(extract_std@a_path, -1)
  expect_lte(extract_std@a_path, 1)

  # Check that they differ from unstandardized
  extract_unstd <- extract_mediation(fit, treatment = "X", mediator = "M")
  expect_false(extract_std@a_path == extract_unstd@a_path)
})

test_that("extract_mediation validates lavaan inputs", {
  skip_if_not_installed("lavaan")

  data <- generate_mediation_data(n = 200, seed = 202)
  model <- "
    M ~ X
    Y ~ M + X
  "
  fit <- lavaan::sem(model, data = data)

  # Invalid treatment
  expect_error(
    extract_mediation(fit, treatment = "NonExistent", mediator = "M"),
    "Treatment variable 'NonExistent' not found"
  )

  # Invalid mediator
  expect_error(
    extract_mediation(fit, treatment = "X", mediator = "NonExistent"),
    "Mediator variable 'NonExistent' not found"
  )

  # Invalid outcome (explicit)
  expect_error(
    extract_mediation(fit, treatment = "X", mediator = "M", outcome = "NonExistent"),
    "Outcome variable 'NonExistent' not found"
  )
})

test_that("extract_mediation auto-detects outcome correctly", {
  skip_if_not_installed("lavaan")

  data <- generate_mediation_data(n = 200, seed = 303)
  model <- "
    M ~ X
    Y ~ M + X
  "
  fit <- lavaan::sem(model, data = data)

  extract <- extract_mediation(fit, treatment = "X", mediator = "M")
  expect_equal(extract@outcome, "Y")

  # Ambiguous outcome case
  # Y1 ~ M + X
  # Y2 ~ M + X
  data$Y2 <- data$Y + rnorm(200)
  model_ambig <- "
    M ~ X
    Y ~ M + X
    Y2 ~ M + X
  "
  fit_ambig <- lavaan::sem(model_ambig, data = data)

  expect_error(
    extract_mediation(fit_ambig, treatment = "X", mediator = "M"),
    "Multiple potential outcomes detected"
  )
})

test_that("lavaan extraction matches lm extraction", {
  skip_if_not_installed("lavaan")

  # Same data, different methods
  data <- generate_mediation_data(n = 500, a = 0.5, b = 0.4, c_prime = 0.2, seed = 404)

  # lavaan
  model <- "
    M ~ X
    Y ~ M + X
  "
  fit_lavaan <- lavaan::sem(model, data = data)
  extract_lavaan <- extract_mediation(fit_lavaan, treatment = "X", mediator = "M")

  # lm
  fit_m <- lm(M ~ X, data = data)
  fit_y <- lm(Y ~ X + M, data = data)
  extract_lm <- extract_mediation(fit_m, model_y = fit_y, treatment = "X", mediator = "M")

  # Should be very similar (ML vs OLS are asymptotically equivalent for Gaussian)
  expect_equal(extract_lavaan@a_path, extract_lm@a_path, tolerance = 0.05)
  expect_equal(extract_lavaan@b_path, extract_lm@b_path, tolerance = 0.05)
  expect_equal(extract_lavaan@c_prime, extract_lm@c_prime, tolerance = 0.05)

  # Residual variances should also match (ML uses n, OLS uses n-p, so small diff)
  expect_equal(extract_lavaan@sigma_m, extract_lm@sigma_m, tolerance = 0.05)
  expect_equal(extract_lavaan@sigma_y, extract_lm@sigma_y, tolerance = 0.05)
})
test_that("extract_mediation handles missing data with FIML", {
  skip_if_not_installed("lavaan")

  data <- generate_mediation_data(n = 200, seed = 505)
  # Introduce missing data
  data$X[1:10] <- NA
  data$M[11:20] <- NA

  model <- "
    M ~ X
    Y ~ M + X
  "

  # Fit with FIML
  fit <- lavaan::sem(model, data = data, missing = "fiml", fixed.x = FALSE)

  extract <- extract_mediation(fit, treatment = "X", mediator = "M")

  # Check n_obs (should be total sample size for FIML, or effective?)
  # lavaan::lavInspect(fit, "nobs") returns number of used observations (patterns)
  # But for FIML it usually returns total N if all rows are used (even partial)
  expect_equal(extract@n_obs, 200)

  # Check data extraction
  # Should return the data used
  expect_equal(nrow(extract@data), 200)
  expect_true(any(is.na(extract@data$X)))
})

test_that("extract_mediation handles multiple mediators", {
  skip_if_not_installed("lavaan")

  data <- generate_mediation_data(n = 200, seed = 606)
  data$M2 <- 0.3 * data$X + rnorm(200)
  data$Y <- 0.4 * data$M + 0.3 * data$M2 + 0.2 * data$X + rnorm(200)

  model <- "
    M ~ a1*X
    M2 ~ a2*X
    Y ~ b1*M + b2*M2 + cp*X
  "

  fit <- lavaan::sem(model, data = data)

  # Extract for M
  extract1 <- extract_mediation(fit, treatment = "X", mediator = "M")
  expect_equal(extract1@mediator, "M")
  expect_equal(extract1@b_path, lavaan::coef(fit)["b1"], ignore_attr = TRUE)

  # Extract for M2
  extract2 <- extract_mediation(fit, treatment = "X", mediator = "M2")
  expect_equal(extract2@mediator, "M2")
  expect_equal(extract2@b_path, lavaan::coef(fit)["b2"], ignore_attr = TRUE)
  expect_equal(extract2@a_path, lavaan::coef(fit)["a2"], ignore_attr = TRUE)
})

test_that("extract_mediation handles robust estimators", {
  skip_if_not_installed("lavaan")

  data <- generate_mediation_data(n = 200, seed = 707)

  model <- "
    M ~ X
    Y ~ M + X
  "

  # Fit with robust estimator
  fit <- lavaan::sem(model, data = data, estimator = "MLR")

  extract <- extract_mediation(fit, treatment = "X", mediator = "M")

  # Check vcov is robust
  vcov_robust <- lavaan::vcov(fit) # This returns robust vcov by default for MLR

  # Verify extracted vcov matches subset of robust vcov
  # We can't easily check "is robust" property on the matrix, but we can check values match
  # the robust vcov from lavaan

  # Re-implement subsetting logic to verify
  pt <- lavaan::parTable(fit)
  a_idx <- pt$free[pt$lhs == "M" & pt$op == "~" & pt$rhs == "X"]

  expect_equal(extract@vcov["a", "a"], vcov_robust[a_idx, a_idx], ignore_attr = TRUE)
})
