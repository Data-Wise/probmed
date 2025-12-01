test_that("pmed works with formula interface", {
  # Simulate data
  set.seed(123)
  n <- 300
  data <- data.frame(
    X = stats::rnorm(n),
    C = stats::rnorm(n)
  )
  data$M <- 0.5 * data$X + 0.3 * data$C + stats::rnorm(n)
  data$Y <- 0.4 * data$M + 0.2 * data$X + 0.2 * data$C + stats::rnorm(n)

  # Compute P_med with plugin
  result <- pmed(
    Y ~ X + M + C,
    formula_m = M ~ X + C,
    data = data,
    treatment = "X",
    mediator = "M",
    method = "plugin"
  )
  # Tests
  expect_true(S7::S7_inherits(result, PmedResult))
  expect_gte(result@estimate, 0)
  expect_lte(result@estimate, 1)
  expect_equal(result@method, "plugin")
})

test_that("parametric bootstrap produces valid CIs", {
  skip_if_not_installed("MASS")

  set.seed(456)
  n <- 200
  data <- data.frame(
    X = stats::rnorm(n),
    C = stats::rnorm(n)
  )
  data$M <- 0.5 * data$X + stats::rnorm(n)
  data$Y <- 0.4 * data$M + 0.2 * data$X + stats::rnorm(n)

  result <- pmed(
    Y ~ X + M + C,
    formula_m = M ~ X + C,
    data = data,
    treatment = "X",
    mediator = "M",
    method = "parametric_bootstrap",
    n_boot = 100, # Small for speed
    seed = 789
  )

  expect_false(is.na(result@ci_lower))
  expect_false(is.na(result@ci_upper))
  expect_lt(result@ci_lower, result@ci_upper)
  expect_length(result@boot_estimates, 100)
})
