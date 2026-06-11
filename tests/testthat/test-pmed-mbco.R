# MBCO method: dispatch, guards, point estimate (CIs added in later tasks)

test_that("method='mbco' dispatches and returns a PmedResult with method label", {
  data <- generate_mediation_data(n = 800, a = 0.5, b = 0.5, c_prime = 0.3, seed = 11)
  res <- pmed(
    Y ~ X + M,
    formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M", method = "mbco"
  )
  expect_s7_class(res, PmedResult)
  expect_identical(res@method, "mbco")
  expect_true(is.na(res@n_boot))
  expect_length(res@boot_estimates, 0)
})

test_that("mbco point estimate matches the plugin/closed-form P_med", {
  data <- generate_mediation_data(n = 4000, a = 0.5, b = 0.5, c_prime = 0.3, seed = 12)
  res <- pmed(Y ~ X + M,
    formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M", method = "mbco"
  )
  plug <- pmed(Y ~ X + M,
    formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M", method = "plugin",
    n_sim = 2e5, seed = 5
  )
  expect_equal(res@estimate, plug@estimate, tolerance = 0.02)
  expect_equal(res@ie_estimate, plug@ie_estimate, tolerance = 1e-8)
})

test_that("mbco errors on a non-Gaussian outcome with an informative message", {
  data <- generate_binary_mediation_data(
    n = 800, a = 0.5, b = 0.8, c_prime = 0.3,
    binary_mediator = FALSE, binary_outcome = TRUE, seed = 21
  )
  expect_error(
    pmed(Y ~ X + M,
      formula_m = M ~ X, data = data,
      treatment = "X", mediator = "M",
      family_y = stats::binomial(), method = "mbco"
    ),
    "Gaussian"
  )
})

test_that("mbco errors when the source data is absent", {
  data <- generate_mediation_data(n = 400, a = 0.5, b = 0.5, c_prime = 0.3, seed = 31)
  fit_m <- stats::glm(M ~ X, data = data)
  fit_y <- stats::glm(Y ~ X + M, data = data)
  ex <- extract_mediation(fit_m,
    treatment = "X", mediator = "M",
    model_y = fit_y, data = data
  )
  ex@data <- NULL
  expect_error(pmed(ex, method = "mbco"), "data")
})

test_that("constrained P_med log-lik at the MLE equals the free log-lik", {
  data <- generate_mediation_data(n = 1500, a = 0.5, b = 0.5, c_prime = 0.3, seed = 41)
  fit_m <- stats::glm(M ~ X, data = data)
  fit_y <- stats::glm(Y ~ X + M, data = data)
  ex <- extract_mediation(fit_m,
    treatment = "X", mediator = "M",
    model_y = fit_y, data = data
  )
  prep <- .mbco_prep(ex, x_ref = 0, x_value = 1)

  # P_med at the free MLE
  a <- ex@a_path
  b <- ex@b_path
  d <- a * b / sqrt(2 * (b^2 * prep$Vm_hat + prep$Vy_hat))
  p_mle <- as.numeric(stats::pnorm(d))

  a_sign <- sign(a)
  if (a_sign == 0) a_sign <- 1
  ll0 <- .mbco_ll_constrained_pmed(prep, qstar = stats::qnorm(p_mle), a_sign = a_sign)

  # Non-binding constraint at the MLE: LR ~ 0, and ll0 never exceeds the free ll.
  expect_lte(ll0, prep$ll_free + 1e-6)
  expect_equal(-2 * (ll0 - prep$ll_free), 0, tolerance = 1e-2)
})

test_that("constrained P_med log-lik handles the null (p*=0.5, b=0) submodel", {
  data <- generate_mediation_data(n = 1000, a = 0.5, b = 0.5, c_prime = 0.3, seed = 42)
  fit_m <- stats::glm(M ~ X, data = data)
  fit_y <- stats::glm(Y ~ X + M, data = data)
  ex <- extract_mediation(fit_m,
    treatment = "X", mediator = "M",
    model_y = fit_y, data = data
  )
  prep <- .mbco_prep(ex, x_ref = 0, x_value = 1)
  ll_null <- .mbco_ll_constrained_pmed(prep, qstar = 0, a_sign = 1)
  # Null nests inside the free model: strictly less likely here (b != 0 truly).
  expect_lt(ll_null, prep$ll_free)
})

test_that("mbco P_med interval brackets the estimate, lies in (0,1), and is deterministic", {
  data <- generate_mediation_data(n = 1200, a = 0.5, b = 0.5, c_prime = 0.3, seed = 51)
  r1 <- pmed(Y ~ X + M,
    formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M", method = "mbco", ci_level = 0.95
  )
  r2 <- pmed(Y ~ X + M,
    formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M", method = "mbco", ci_level = 0.95
  )
  expect_gt(r1@ci_lower, 0)
  expect_lt(r1@ci_upper, 1)
  expect_lte(r1@ci_lower, r1@estimate)
  expect_gte(r1@ci_upper, r1@estimate)
  expect_equal(r1@ci_level, 0.95)
  # Seed-free: identical across calls.
  expect_identical(r1@ci_lower, r2@ci_lower)
  expect_identical(r1@ci_upper, r2@ci_upper)
})

test_that("profiled mbco interval matches the independent oracle (no covariates)", {
  data <- generate_mediation_data(n = 1500, a = 0.5, b = 0.5, c_prime = 0.3, seed = 61)
  res <- pmed(Y ~ X + M,
    formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M", method = "mbco"
  )
  orc <- oracle_mbco_ci(data, level = 0.95)
  expect_equal(res@estimate, unname(orc["estimate"]), tolerance = 1e-6)
  expect_equal(res@ci_lower, unname(orc["lower"]), tolerance = 0.01)
  expect_equal(res@ci_upper, unname(orc["upper"]), tolerance = 0.01)
})

test_that("profiled mbco interval matches the independent oracle (one covariate)", {
  data <- generate_mediation_data(
    n = 1500, a = 0.5, b = 0.5, c_prime = 0.3,
    n_covariates = 1, covariate_effects = list(m = 0.3, y = 0.4),
    seed = 62
  )
  res <- pmed(Y ~ X + M + C1,
    formula_m = M ~ X + C1, data = data,
    treatment = "X", mediator = "M", method = "mbco"
  )
  orc <- oracle_mbco_ci(data, level = 0.95, covs = "C1")
  expect_equal(res@estimate, unname(orc["estimate"]), tolerance = 1e-6)
  expect_equal(res@ci_lower, unname(orc["lower"]), tolerance = 0.015)
  expect_equal(res@ci_upper, unname(orc["upper"]), tolerance = 0.015)
})

test_that("mbco IE interval brackets a*b and excludes 0 for a strong effect", {
  data <- generate_mediation_data(n = 1500, a = 0.6, b = 0.6, c_prime = 0.3, seed = 71)
  res <- pmed(Y ~ X + M,
    formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M", method = "mbco"
  )
  expect_lte(res@ie_ci_lower, res@ie_estimate)
  expect_gte(res@ie_ci_upper, res@ie_estimate)
  expect_gt(res@ie_ci_lower, 0) # strong positive mediation
  expect_lt(res@ie_ci_lower, res@ie_ci_upper)
})

test_that("mbco IE interval includes 0 at a near-null indirect effect", {
  data <- generate_mediation_data(n = 1500, a = 0.02, b = 0.02, c_prime = 0.3, seed = 72)
  res <- pmed(Y ~ X + M,
    formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M", method = "mbco"
  )
  expect_lt(res@ie_ci_lower, 0)
  expect_gt(res@ie_ci_upper, 0)
})

test_that("mbco IE constrained log-lik finds the global optimum on the bimodal tail", {
  data <- generate_mediation_data(n = 1500, a = 0.6, b = 0.6, c_prime = 0.3, seed = 71)
  fit_m <- stats::glm(M ~ X, data = data)
  fit_y <- stats::glm(Y ~ X + M, data = data)
  ex <- extract_mediation(fit_m,
    treatment = "X", mediator = "M",
    model_y = fit_y, data = data
  )
  prep <- .mbco_prep(ex, x_ref = 0, x_value = 1)
  ie0 <- 0.05 # deep in the tail, where the a-profile is bimodal

  # Independent brute-force global maximum: profile over a on a dense grid,
  # with Vm, Vy, and all linear nuisance terms in closed form (MLE = sse/n).
  n <- prep$n
  grid_ll <- function(a) {
    if (abs(a) < 1e-6) {
      return(-Inf)
    }
    b <- ie0 / a
    sse_m <- sum(stats::lm.fit(prep$Dm, prep$M - a * prep$X)$residuals^2)
    sse_y <- sum(stats::lm.fit(prep$Dy, prep$Y - b * prep$M)$residuals^2)
    -n / 2 * log(2 * pi * sse_m / n) - n / 2 -
      n / 2 * log(2 * pi * sse_y / n) - n / 2
  }
  a_grid <- seq(0.01, 3.0, length.out = 6000)
  brute <- max(vapply(a_grid, grid_ll, numeric(1)))

  est <- .mbco_ll_constrained_ie(prep, ie0)
  expect_equal(est, brute, tolerance = 1e-3)
})

test_that("mbco P_med interval reflects correctly under a reversed treatment contrast", {
  # Reversing the contrast (delta -> -delta) maps P_med -> 1 - P_med, so the
  # whole interval should reflect: [L, U] at delta=+1 becomes [1-U, 1-L] at delta=-1.
  data <- generate_mediation_data(n = 1200, a = 0.5, b = 0.5, c_prime = 0.3, seed = 51)
  fwd <- pmed(Y ~ X + M,
    formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M", method = "mbco",
    x_ref = 0, x_value = 1
  )
  rev <- pmed(Y ~ X + M,
    formula_m = M ~ X, data = data,
    treatment = "X", mediator = "M", method = "mbco",
    x_ref = 1, x_value = 0
  )
  expect_equal(rev@estimate, 1 - fwd@estimate, tolerance = 1e-6)
  expect_equal(rev@ci_lower, 1 - fwd@ci_upper, tolerance = 1e-3)
  expect_equal(rev@ci_upper, 1 - fwd@ci_lower, tolerance = 1e-3)
})

test_that("mbco errors when x_ref equals x_value (degenerate contrast)", {
  data <- generate_mediation_data(n = 400, a = 0.5, b = 0.5, c_prime = 0.3, seed = 31)
  expect_error(
    pmed(Y ~ X + M,
      formula_m = M ~ X, data = data,
      treatment = "X", mediator = "M", method = "mbco",
      x_ref = 1, x_value = 1
    ),
    "x_value|x_ref|contrast"
  )
})

test_that(".mbco_invert pins to the domain bounds when the statistic never crosses", {
  flat <- function(t) -1 # excess always negative: never rejected
  ci <- .mbco_invert(0.5, flat, domain = c(0.01, 0.99), step = 0.01)
  expect_equal(unname(ci["lower"]), 0.01)
  expect_equal(unname(ci["upper"]), 0.99)
})

test_that(".mbco_invert returns NA when the estimate is already rejected", {
  always_pos <- function(t) 1 # excess > 0 immediately: same-sign bracket
  ci <- .mbco_invert(0.5, always_pos, domain = c(0.01, 0.99), step = 0.01)
  expect_true(is.na(ci["lower"]))
  expect_true(is.na(ci["upper"]))
})
