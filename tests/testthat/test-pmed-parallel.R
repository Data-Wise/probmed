# Parallel (joint) P_med: dispatch, estimand, four methods, guards

# ---- Dispatch & basic structure ----

test_that("pmed(ParallelMediationData) dispatches and labels the method", {
  data <- generate_parallel_data(n = 800, seed = 101)
  ex <- extract_parallel(data)
  res <- pmed(ex, method = "plugin")
  expect_s7_class(res, PmedResult)
  expect_identical(res@method, "plugin")
  expect_true(is.na(res@n_boot))
  expect_length(res@boot_estimates, 0)
})

# ---- Estimand correctness ----

test_that(".pmed_parallel_closed reduces to the single-mediator formula at k=1", {
  a <- 0.5; b <- 0.6; sm <- 0.9; sy <- 1.1; delta <- 1
  single <- as.numeric(stats::pnorm(a * b / sqrt(2 * (b^2 * sm^2 + sy^2))))
  par <- .pmed_parallel_closed(a, b, sm, sy, delta)
  expect_equal(par, single, tolerance = 1e-12)
})

test_that("plugin joint P_med matches the independent brute-force simulator", {
  data <- generate_parallel_data(
    n = 4000, a_vec = c(0.5, 0.4), b_vec = c(0.5, 0.6),
    c_prime = 0.3, seed = 102
  )
  ex <- extract_parallel(data)
  res <- pmed(ex, method = "plugin")
  orc <- oracle_parallel_pmed(
    ex@a_paths, ex@b_paths, ex@sigma_mediators, ex@sigma_y,
    delta = 1, n_sim = 3e6, seed = 7
  )
  expect_equal(res@estimate, orc, tolerance = 0.01)
})

test_that("joint P_med is invariant to the direct effect c'", {
  d0 <- generate_parallel_data(n = 1500, c_prime = 0.0, seed = 103)
  d5 <- generate_parallel_data(n = 1500, c_prime = 5.0, seed = 103)
  e0 <- pmed(extract_parallel(d0), method = "plugin")@estimate
  e5 <- pmed(extract_parallel(d5), method = "plugin")@estimate
  expect_equal(e0, e5, tolerance = 1e-8)
})

test_that("indirect effect equals sum(a_j b_j)", {
  data <- generate_parallel_data(n = 1200, seed = 104)
  ex <- extract_parallel(data)
  res <- pmed(ex, method = "plugin")
  expect_equal(res@ie_estimate, sum(ex@a_paths * ex@b_paths), tolerance = 1e-10)
})

test_that("reversing the contrast reflects joint P_med to 1 - P_med", {
  data <- generate_parallel_data(n = 1500, seed = 105)
  ex <- extract_parallel(data)
  fwd <- pmed(ex, method = "plugin", x_ref = 0, x_value = 1)
  rev <- pmed(ex, method = "plugin", x_ref = 1, x_value = 0)
  expect_equal(rev@estimate, 1 - fwd@estimate, tolerance = 1e-10)
})

test_that("three parallel mediators agree with the simulator", {
  data <- generate_parallel_data(
    n = 4000, a_vec = c(0.4, 0.3, 0.5), b_vec = c(0.5, 0.4, 0.3),
    seed = 106
  )
  ex <- extract_parallel(data)
  res <- pmed(ex, method = "plugin")
  orc <- oracle_parallel_pmed(
    ex@a_paths, ex@b_paths, ex@sigma_mediators, ex@sigma_y,
    delta = 1, n_sim = 3e6, seed = 9
  )
  expect_equal(res@estimate, orc, tolerance = 0.012)
})

# ---- Guards ----

test_that("non-Gaussian (missing residual SD) parallel extract errors", {
  data <- generate_parallel_data(n = 600, seed = 107)
  ex <- extract_parallel(data)
  ex@sigma_y <- numeric(0) # mark "not supplied" (non-Gaussian indicator)
  expect_error(pmed(ex, method = "plugin"), "Gaussian")
})

test_that("nonparametric bootstrap errors when data is absent", {
  data <- generate_parallel_data(n = 600, seed = 108)
  ex <- extract_parallel(data)
  ex@data <- NULL
  expect_error(pmed(ex, method = "nonparametric_bootstrap"), "[Dd]ata")
})

# ---- Parametric bootstrap ----

test_that("parametric bootstrap CI brackets the estimate, lies in (0,1), is seeded", {
  data <- generate_parallel_data(n = 1200, seed = 109)
  ex <- extract_parallel(data)
  r1 <- pmed(ex, method = "parametric_bootstrap", n_boot = 600, seed = 1)
  r2 <- pmed(ex, method = "parametric_bootstrap", n_boot = 600, seed = 1)
  expect_gt(r1@ci_lower, 0)
  expect_lt(r1@ci_upper, 1)
  expect_lte(r1@ci_lower, r1@estimate)
  expect_gte(r1@ci_upper, r1@estimate)
  expect_lt(r1@ci_lower, r1@ci_upper)
  expect_identical(r1@ci_lower, r2@ci_lower)
  expect_identical(r1@ci_upper, r2@ci_upper)
})

test_that("parametric bootstrap IE interval excludes 0 for a strong joint effect", {
  data <- generate_parallel_data(
    n = 1500, a_vec = c(0.6, 0.6), b_vec = c(0.6, 0.6), seed = 110
  )
  ex <- extract_parallel(data)
  res <- pmed(ex, method = "parametric_bootstrap", n_boot = 800, seed = 2)
  expect_gt(res@ie_ci_lower, 0)
  expect_lte(res@ie_ci_lower, res@ie_estimate)
  expect_gte(res@ie_ci_upper, res@ie_estimate)
})

# ---- Nonparametric bootstrap ----

test_that("nonparametric bootstrap CI brackets the estimate and is seeded", {
  data <- generate_parallel_data(n = 1200, seed = 111)
  ex <- extract_parallel(data)
  r1 <- pmed(ex, method = "nonparametric_bootstrap", n_boot = 400, seed = 3)
  r2 <- pmed(ex, method = "nonparametric_bootstrap", n_boot = 400, seed = 3)
  expect_gt(r1@ci_lower, 0)
  expect_lt(r1@ci_upper, 1)
  expect_lte(r1@ci_lower, r1@estimate)
  expect_gte(r1@ci_upper, r1@estimate)
  expect_identical(r1@ci_lower, r2@ci_lower)
})

# ---- MBCO ----

test_that("mbco point estimate matches the plugin/closed-form joint P_med", {
  data <- generate_parallel_data(n = 2000, seed = 112)
  ex <- extract_parallel(data)
  m <- pmed(ex, method = "mbco")
  p <- pmed(ex, method = "plugin")
  expect_equal(m@estimate, p@estimate, tolerance = 1e-8)
  expect_equal(m@ie_estimate, p@ie_estimate, tolerance = 1e-8)
})

test_that("mbco joint interval brackets the estimate, lies in (0,1), deterministic", {
  data <- generate_parallel_data(n = 1500, seed = 113)
  ex <- extract_parallel(data)
  r1 <- pmed(ex, method = "mbco")
  r2 <- pmed(ex, method = "mbco")
  expect_gt(r1@ci_lower, 0)
  expect_lt(r1@ci_upper, 1)
  expect_lte(r1@ci_lower, r1@estimate)
  expect_gte(r1@ci_upper, r1@estimate)
  expect_identical(r1@ci_lower, r2@ci_lower)
  expect_identical(r1@ci_upper, r2@ci_upper)
})

test_that("constrained joint P_med log-lik at the MLE equals the free log-lik", {
  data <- generate_parallel_data(n = 1500, seed = 114)
  ex <- extract_parallel(data)
  prep <- .mbco_prep_parallel(ex, x_ref = 0, x_value = 1)
  s_hat <- sum(prep$a_hat * prep$b_hat)
  d <- prep$delta * s_hat /
    sqrt(2 * sum(prep$b_hat^2 * prep$Vm_hat) + 2 * prep$Vy_hat)
  p_mle <- as.numeric(stats::pnorm(d))
  ll0 <- .mbco_ll_constrained_pmed_parallel(prep, stats::qnorm(p_mle))
  expect_lte(ll0, prep$ll_free + 1e-6)
  expect_equal(-2 * (ll0 - prep$ll_free), 0, tolerance = 1e-2)
})

test_that("profiled mbco joint interval matches the full-dimensional oracle", {
  data <- generate_parallel_data(n = 1500, seed = 115)
  ex <- extract_parallel(data)
  res <- pmed(ex, method = "mbco")
  orc <- oracle_parallel_mbco_ci(data, meds = c("M1", "M2"), delta = 1)
  # The estimate assertion is a sanity guard (both paths share the closed form by
  # construction); the CI endpoints are the substantive check -- they come from
  # an independent un-profiled full-dimensional LR inversion in the oracle.
  expect_equal(res@estimate, unname(orc["estimate"]), tolerance = 1e-6)
  expect_equal(res@ci_lower, unname(orc["lower"]), tolerance = 0.02)
  expect_equal(res@ci_upper, unname(orc["upper"]), tolerance = 0.02)
})

test_that("mbco joint interval reflects under a reversed contrast", {
  data <- generate_parallel_data(n = 1500, seed = 113)
  ex <- extract_parallel(data)
  fwd <- pmed(ex, method = "mbco", x_ref = 0, x_value = 1)
  rev <- pmed(ex, method = "mbco", x_ref = 1, x_value = 0)
  expect_equal(rev@estimate, 1 - fwd@estimate, tolerance = 1e-6)
  expect_equal(rev@ci_lower, 1 - fwd@ci_upper, tolerance = 1e-2)
  expect_equal(rev@ci_upper, 1 - fwd@ci_lower, tolerance = 1e-2)
})

test_that("mbco IE interval brackets sum(a_j b_j) and excludes 0 for a strong effect", {
  data <- generate_parallel_data(
    n = 1500, a_vec = c(0.6, 0.5), b_vec = c(0.6, 0.5), seed = 116
  )
  ex <- extract_parallel(data)
  res <- pmed(ex, method = "mbco")
  expect_lte(res@ie_ci_lower, res@ie_estimate)
  expect_gte(res@ie_ci_upper, res@ie_estimate)
  expect_gt(res@ie_ci_lower, 0)
  expect_lt(res@ie_ci_lower, res@ie_ci_upper)
})

test_that("mbco errors on a degenerate contrast (x_ref == x_value)", {
  data <- generate_parallel_data(n = 600, seed = 117)
  ex <- extract_parallel(data)
  expect_error(
    pmed(ex, method = "mbco", x_ref = 1, x_value = 1),
    "x_value|x_ref|contrast"
  )
})

test_that("mbco S=0 (p*=0.5) fit is continuous with the qstar->0 limit (k>=2)", {
  # sum(a*b) ~ 0 with all b_j != 0 (the mediators cancel). The constrained fit at
  # p* = 0.5 must allow b != 0 on the surface S = 0, not force b = 0 -- the latter
  # was ~80-160 log-lik units worse, a discontinuity that spuriously inflated the
  # LR statistic near P_med = 0.5.
  data <- generate_parallel_data(
    n = 1500, a_vec = c(0.4, 0.4), b_vec = c(0.35, -0.30), seed = 121
  )
  ex <- extract_parallel(data)
  prep <- .mbco_prep_parallel(ex, 0, 1)
  ll_at <- .mbco_ll_constrained_pmed_parallel(prep, stats::qnorm(0.5))
  ll_near <- .mbco_ll_constrained_pmed_parallel(prep, stats::qnorm(0.5 + 1e-4))
  expect_equal(ll_at, ll_near, tolerance = 1.0)
})

test_that("profiled mbco joint interval matches the oracle for three mediators", {
  data <- generate_parallel_data(
    n = 1000, a_vec = c(0.4, 0.3, 0.5), b_vec = c(0.5, 0.4, 0.3), seed = 118
  )
  ex <- extract_parallel(data)
  res <- pmed(ex, method = "mbco")
  orc <- oracle_parallel_mbco_ci(data, meds = c("M1", "M2", "M3"), delta = 1)
  expect_equal(res@estimate, unname(orc["estimate"]), tolerance = 1e-6)
  expect_equal(res@ci_lower, unname(orc["lower"]), tolerance = 0.025)
  expect_equal(res@ci_upper, unname(orc["upper"]), tolerance = 0.025)
})

test_that("all methods handle covariates in the parallel mediator/outcome models", {
  data <- generate_parallel_data(n = 1500, n_cov = 1, seed = 119)
  ex <- extract_parallel(data, covariates = "C1")
  plug <- pmed(ex, method = "plugin")
  par <- pmed(ex, method = "parametric_bootstrap", n_boot = 400, seed = 1)
  npar <- pmed(ex, method = "nonparametric_bootstrap", n_boot = 200, seed = 1)
  mb <- pmed(ex, method = "mbco")
  # Plugin point estimate matches the simulator built from the covariate-adjusted
  # path estimates.
  orc <- oracle_parallel_pmed(
    ex@a_paths, ex@b_paths, ex@sigma_mediators, ex@sigma_y,
    delta = 1, n_sim = 3e6, seed = 8
  )
  expect_equal(plug@estimate, orc, tolerance = 0.012)
  expect_equal(par@estimate, plug@estimate, tolerance = 0.03)
  expect_equal(npar@estimate, plug@estimate, tolerance = 0.03)
  # MBCO interval valid and brackets the estimate with covariates present.
  expect_gt(mb@ci_lower, 0)
  expect_lt(mb@ci_upper, 1)
  expect_lte(mb@ci_lower, mb@estimate)
  expect_gte(mb@ci_upper, mb@estimate)
})

test_that("parametric bootstrap works on a lavaan-sourced parallel extract", {
  skip_if_not_installed("lavaan")
  data <- generate_parallel_data(n = 1200, seed = 120)
  ex <- extract_parallel_lavaan(data, meds = c("M1", "M2"))
  skip_if(is.null(ex), "lavaan parallel extraction unavailable")
  # The structural aliases a{j}/b{j} (present in both lm and lavaan extracts)
  # make the parametric draw source-agnostic; this errored before the fix.
  res <- pmed(ex, method = "parametric_bootstrap", n_boot = 400, seed = 1)
  expect_s7_class(res, PmedResult)
  expect_gt(res@ci_lower, 0)
  expect_lt(res@ci_upper, 1)
  expect_lte(res@ci_lower, res@ci_upper)
})

test_that("Gaussian guard fires on a non-Gaussian mediator and on non-plugin entry", {
  data <- generate_parallel_data(n = 600, seed = 122)
  ex <- extract_parallel(data)
  ex@sigma_mediators <- c(ex@sigma_mediators[1], NA_real_) # one non-Gaussian M
  expect_error(pmed(ex, method = "plugin"), "Gaussian")
  expect_error(pmed(ex, method = "parametric_bootstrap"), "Gaussian")
  expect_error(pmed(ex, method = "mbco"), "Gaussian")
})
