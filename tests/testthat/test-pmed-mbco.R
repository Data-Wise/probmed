# MBCO method: dispatch, guards, point estimate (CIs added in later tasks)

test_that("method='mbco' dispatches and returns a PmedResult with method label", {
  data <- generate_mediation_data(n = 800, a = 0.5, b = 0.5, c_prime = 0.3, seed = 11)
  res <- pmed(
    Y ~ X + M, formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M", method = "mbco"
  )
  expect_s7_class(res, PmedResult)
  expect_identical(res@method, "mbco")
  expect_true(is.na(res@n_boot))
  expect_length(res@boot_estimates, 0)
})

test_that("mbco point estimate matches the plugin/closed-form P_med", {
  data <- generate_mediation_data(n = 4000, a = 0.5, b = 0.5, c_prime = 0.3, seed = 12)
  res  <- pmed(Y ~ X + M, formula_m = M ~ X, data = data,
               treatment = "X", mediator = "M", method = "mbco")
  plug <- pmed(Y ~ X + M, formula_m = M ~ X, data = data,
               treatment = "X", mediator = "M", method = "plugin",
               n_sim = 2e5, seed = 5)
  expect_equal(res@estimate, plug@estimate, tolerance = 0.02)
  expect_equal(res@ie_estimate, plug@ie_estimate, tolerance = 1e-8)
})

test_that("mbco errors on a non-Gaussian outcome with an informative message", {
  data <- generate_binary_mediation_data(
    n = 800, a = 0.5, b = 0.8, c_prime = 0.3,
    binary_mediator = FALSE, binary_outcome = TRUE, seed = 21
  )
  expect_error(
    pmed(Y ~ X + M, formula_m = M ~ X, data = data,
         treatment = "X", mediator = "M",
         family_y = stats::binomial(), method = "mbco"),
    "Gaussian"
  )
})

test_that("mbco errors when the source data is absent", {
  data <- generate_mediation_data(n = 400, a = 0.5, b = 0.5, c_prime = 0.3, seed = 31)
  fit_m <- stats::glm(M ~ X, data = data)
  fit_y <- stats::glm(Y ~ X + M, data = data)
  ex <- extract_mediation(fit_m, treatment = "X", mediator = "M",
                          model_y = fit_y, data = data)
  ex@data <- NULL
  expect_error(pmed(ex, method = "mbco"), "data")
})
