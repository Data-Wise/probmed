# Tests pinning P_med to the *mediation* estimand (Definition 1):
#   P_med = P(Y(x, M(x)) > Y(x, M(x*))) + 0.5 P(Y(x, M(x)) = Y(x, M(x*)))
# Treatment is held at x for both outcomes (so c' cancels); the mediator is
# drawn under both levels independently. See planning/SPEC-estimand-fix.md.

# Closed form for a Gaussian DGM: Phi(ab / sqrt(2 (b^2 sigma_m^2 + sigma_y^2)))
pmed_closed_form <- function(a, b, sigma_m = 1, sigma_y = 1) {
  stats::pnorm(a * b / sqrt(2 * (b^2 * sigma_m^2 + sigma_y^2)))
}

# ==============================================================================
# Core estimand: closed-form agreement + c' invariance (deterministic)
# ==============================================================================

test_that(".pmed_core_simple matches the Gaussian closed form", {
  a <- 0.5; b <- 0.5; sm <- 1; sy <- 1
  target <- pmed_closed_form(a, b, sm, sy) # ~0.5628
  set.seed(101)
  est <- .pmed_core_simple(
    a = a, b = b, c_prime = 0.3, x_ref = 0, x_value = 1,
    sigma_y = sy, sigma_m = sm, n_sim = 2e5
  )
  expect_equal(est, target, tolerance = 0.01)
  expect_gt(est, 0.5) # positive mediation when a, b > 0
})

test_that("P_med is invariant to the direct effect c' (Gaussian)", {
  # Same RNG stream + same a, b, sigma => identical draws; only c' differs.
  # Because both outcomes hold X = x_value, c' cancels exactly.
  core <- function(cp) {
    set.seed(202)
    .pmed_core_simple(
      a = 0.5, b = 0.5, c_prime = cp, x_ref = 0, x_value = 1,
      sigma_y = 1, sigma_m = 1, n_sim = 5e4
    )
  }
  expect_equal(core(0.0), core(0.6))
  expect_equal(core(-0.4), core(0.9))
})

test_that("the old direct-effect contrast is NOT what we compute", {
  # Guard against regression to P(Y(x*,M_x) > Y(x,M_x)) ~ Phi(-c'/sqrt2) = 0.416.
  set.seed(303)
  est <- .pmed_core_simple(
    a = 0.5, b = 0.5, c_prime = 0.3, x_ref = 0, x_value = 1,
    sigma_y = 1, sigma_m = 1, n_sim = 2e5
  )
  expect_gt(est, 0.5) # mediation value is > 0.5, not the 0.416 direct contrast
})

# ==============================================================================
# Integration through the formula interface (Gaussian + binary)
# ==============================================================================

test_that("plugin P_med (Gaussian GLM) reproduces the closed form", {
  data <- generate_mediation_data(n = 4000, a = 0.5, b = 0.5, c_prime = 0.3, seed = 11)
  res <- pmed(
    Y ~ X + M, formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M", method = "plugin", n_sim = 2e5, seed = 5
  )
  expect_equal(res@estimate, pmed_closed_form(0.5, 0.5), tolerance = 0.02)
})

test_that("binary-outcome P_med is non-degenerate and exceeds 0.5", {
  data <- generate_binary_mediation_data(
    n = 4000, a = 0.5, b = 0.8, c_prime = 0.3,
    binary_mediator = FALSE, binary_outcome = TRUE, seed = 21
  )
  res <- pmed(
    Y ~ X + M, formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M", family_y = stats::binomial(),
    method = "plugin", n_sim = 2e5, seed = 5
  )
  expect_gt(res@estimate, 0) # strictly interior, not the degenerate 0/1
  expect_lt(res@estimate, 1)
  expect_gt(res@estimate, 0.5)
})

# ==============================================================================
# Bootstraps: IE must equal a*b (the parametric indexing was previously wrong)
# ==============================================================================

test_that("parametric bootstrap IE equals a*b and P_med is sane", {
  skip_if_not_installed("MASS")
  data <- generate_mediation_data(n = 1000, a = 0.5, b = 0.4, c_prime = 0.3, seed = 31)
  res <- pmed(
    Y ~ X + M, formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M",
    method = "parametric_bootstrap", n_boot = 400, seed = 7
  )
  plug <- pmed(Y ~ X + M, formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M", method = "plugin", n_sim = 2e5, seed = 7)
  # IE bootstrap mean must track a*b (not the intercept product it used to)
  expect_equal(res@ie_estimate, plug@ie_estimate, tolerance = 0.05)
  expect_gt(res@estimate, 0.5)
  expect_lt(res@ci_lower, res@ci_upper)
})

test_that("nonparametric bootstrap IE equals the fitted a*b", {
  data <- generate_mediation_data(n = 1000, a = 0.5, b = 0.4, c_prime = 0.3, seed = 41)
  fitted_ie <- unname(coef(glm(M ~ X, data = data))["X"] *
    coef(glm(Y ~ X + M, data = data))["M"])
  res <- pmed(
    Y ~ X + M, formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M",
    method = "nonparametric_bootstrap", n_boot = 300, seed = 7
  )
  # bootstrap-mean IE must track the fitted a*b (resampling noise only)
  expect_equal(res@ie_estimate, fitted_ie, tolerance = 0.1)
  expect_gt(res@estimate, 0.5)
})
