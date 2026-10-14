# Opt-in model-based response bootstrap (repair tasks T12, T13).

mb_gauss <- function(n = 300, seed = 41) {
  set.seed(seed)
  X <- stats::rbinom(n, 1, 0.5)
  C <- stats::rnorm(n)
  M <- 0.5 * X + 0.4 * C + stats::rnorm(n)
  Y <- 0.3 * X + 0.4 * M + 0.4 * C + stats::rnorm(n)
  data.frame(X = X, C = C, M = M, Y = Y)
}

# Independent reference: delta-method SE of the closed-form Gaussian P_med
# g = Phi(a b / sqrt(2 (b^2 Vm + Vy))), with theta = (a, b, Vm, Vy). Under the
# Gaussian sequential model a_hat, b_hat, Vm_hat, Vy_hat are mutually
# independent; Var(V_hat) = 2 V^2 / df.
delta_se_pmed <- function(d) {
  fm <- stats::lm(M ~ X + C, d)
  fy <- stats::lm(Y ~ X + M + C, d)
  a <- stats::coef(fm)[["X"]]; b <- stats::coef(fy)[["M"]]
  Vm <- summary(fm)$sigma^2; Vy <- summary(fy)$sigma^2
  g <- function(t) stats::pnorm(t[1] * t[2] / sqrt(2 * (t[2]^2 * t[3] + t[4])))
  th <- c(a, b, Vm, Vy)
  grad <- vapply(1:4, function(j) {
    h <- 1e-6 * max(1, abs(th[j])); e <- replace(numeric(4), j, h)
    (g(th + e) - g(th - e)) / (2 * h)
  }, numeric(1))
  v <- c(stats::vcov(fm)["X", "X"], stats::vcov(fy)["M", "M"],
         2 * Vm^2 / fm$df.residual, 2 * Vy^2 / fy$df.residual)
  sqrt(sum(grad^2 * v))
}

run_mb <- function(d, ...) {
  pmed(Y ~ X + M + C, formula_m = M ~ X + C, data = d, treatment = "X",
       mediator = "M", method = "model_bootstrap", ...)
}

test_that("bootstrap SD agrees with the delta-method SE (known answer)", {
  skip_on_cran()
  d <- mb_gauss()
  ref <- delta_se_pmed(d)
  r <- run_mb(d, n_boot = 400, n_sim = 20000, seed = 1)
  boot_sd <- stats::sd(r@boot_estimates)
  expect_lt(abs(boot_sd / ref - 1), 0.20)
  expect_identical(r@diagnostics$interval_status, "ok")
  expect_identical(unname(r@diagnostics$endpoint_status), c("percentile", "percentile"))
  # Positive control: draws without refitting (plug-in noise only) are far
  # narrower than the reference, so the 20% band can fail.
  noise_sd <- stats::sd(vapply(1:40, function(s) pmed(
    Y ~ X + M + C, formula_m = M ~ X + C, data = d, treatment = "X",
    mediator = "M", method = "plugin", n_sim = 20000, seed = s)@estimate,
    numeric(1)))
  expect_gt(abs(noise_sd / ref - 1), 0.5)
})

test_that("the point estimate is the fitted-model plug-in, not the draw mean", {
  d <- mb_gauss()
  r <- run_mb(d, n_boot = 30, n_sim = 5000, seed = 9)
  p <- pmed(Y ~ X + M + C, formula_m = M ~ X + C, data = d, treatment = "X",
            mediator = "M", method = "plugin", n_sim = 500000, seed = 9)
  # The point estimate uses 100 * n_sim pairs (T21: its MC error is not
  # averaged away by the percentile interval).
  expect_identical(r@estimate, p@estimate)
  expect_false(isTRUE(all.equal(r@estimate, mean(r@boot_estimates))))
})

test_that("fewer than 99% valid draws give an incomplete interval, never recycled draws", {
  d <- mb_gauss(n = 150)
  core <- .pmed_core_simple
  calls <- 0L
  fail_every <- function(k) function(...) {
    calls <<- calls + 1L
    # call 1 is the point estimate; fail every k-th draw after it
    if (calls > 1L && (calls - 1L) %% k == 0L) NA_real_ else core(...)
  }
  testthat::local_mocked_bindings(.pmed_core_simple = fail_every(50))
  calls <- 0L
  r <- run_mb(d, n_boot = 100, n_sim = 1000, seed = 2)   # 2 of 100 fail
  expect_identical(r@diagnostics$interval_status, "bootstrap_incomplete")
  expect_true(is.na(r@ci_lower) && is.na(r@ci_upper))
  expect_identical(r@diagnostics$n_valid, 98L)
  expect_identical(r@diagnostics$failures[["nonfinite"]], 2L)
  expect_identical(sum(is.na(r@boot_estimates)), 2L)
  expect_false(r@converged)
})

test_that("exactly 99% valid draws still report the interval", {
  d <- mb_gauss(n = 150)
  core <- .pmed_core_simple
  calls <- 0L
  testthat::local_mocked_bindings(.pmed_core_simple = function(...) {
    calls <<- calls + 1L
    if (calls == 51L) NA_real_ else core(...)
  })
  r <- run_mb(d, n_boot = 100, n_sim = 1000, seed = 2)
  expect_identical(r@diagnostics$n_valid, 99L)
  expect_identical(r@diagnostics$interval_status, "ok")
  expect_true(r@ci_lower < r@ci_upper)
})

test_that("separated draws are counted and reported as failures", {
  skip_if_not_installed("lpSolve")
  set.seed(12)
  n <- 40
  X <- stats::rbinom(n, 1, 0.5)
  M <- 0.5 * X + stats::rnorm(n)
  Y <- stats::rbinom(n, 1, stats::plogis(-2.2 + 0.8 * M))
  d <- data.frame(X = X, M = M, Y = Y)
  r <- suppressWarnings(pmed(Y ~ X + M, formula_m = M ~ X, data = d,
                             treatment = "X", mediator = "M",
                             family_y = stats::binomial(),
                             method = "model_bootstrap", n_boot = 60,
                             n_sim = 2000, seed = 3))
  expect_gt(r@diagnostics$failures[["separated"]], 0L)
  expect_identical(r@diagnostics$interval_status, "bootstrap_incomplete")
})

test_that("binary models with covariates are standardized in every draw", {
  set.seed(13)
  n <- 300
  X <- stats::rbinom(n, 1, 0.5)
  C <- stats::rnorm(n)
  M <- stats::rbinom(n, 1, stats::plogis(-0.25 + 0.8 * X + 0.6 * C))
  Y <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.3 * X + 1.0 * M + 0.8 * C))
  d <- data.frame(X = X, C = C, M = M, Y = Y)
  r <- pmed(Y ~ X + M + C, formula_m = M ~ X + C, data = d, treatment = "X",
            mediator = "M", family_m = stats::binomial(),
            family_y = stats::binomial(), method = "model_bootstrap",
            n_boot = 50, n_sim = 4000, seed = 4)
  expect_identical(r@diagnostics$interval_status, "ok")
  expect_true(r@ci_lower <= r@estimate && r@estimate <= r@ci_upper)
})

test_that("model_bootstrap is opt-in: the default method is unchanged", {
  m <- S7::method(pmed, S7::class_formula)
  expect_identical(eval(formals(m)$method)[1], "parametric_bootstrap")
})

test_that("each draw-failure reason is counted, never hidden", {
  d <- mb_gauss(n = 120)
  # Fit outside the mocks: pmed()'s own glm() calls also use glm.fit.
  ex <- extract_mediation(stats::glm(M ~ X + C, data = d), treatment = "X",
                          mediator = "M", model_y = stats::glm(Y ~ X + M + C, data = d),
                          data = d)
  run <- function() .pmed_model_boot(ex, 0, 1, n_boot = 10, ci_level = 0.95,
                                     seed = 5, n_sim = 500)
  real_fit <- stats::glm.fit
  # Fit error.
  local({
    testthat::local_mocked_bindings(glm.fit = function(...) stop("boom"), .package = "stats")
    r <- run()
    expect_identical(r@diagnostics$failures[["fit_error"]], 10L)
    expect_identical(r@diagnostics$interval_status, "bootstrap_incomplete")
  })
  # Non-convergence.
  local({
    testthat::local_mocked_bindings(glm.fit = function(...) {
      f <- real_fit(...); f$converged <- FALSE; f
    }, .package = "stats")
    expect_identical(run()@diagnostics$failures[["nonconvergence"]], 10L)
  })
  # Non-finite coefficients.
  local({
    testthat::local_mocked_bindings(glm.fit = function(...) {
      f <- real_fit(...); f$coefficients[1] <- NA_real_; f
    }, .package = "stats")
    expect_identical(run()@diagnostics$failures[["nonfinite"]], 10L)
  })
})

test_that("separated mediator draws are counted, and source data is required", {
  set.seed(14)
  n <- 150
  X <- stats::rbinom(n, 1, 0.5)
  M <- stats::rbinom(n, 1, stats::plogis(-0.2 + 0.8 * X))
  d <- data.frame(X = X, M = M, Y = 0.5 * M + stats::rnorm(n))
  ex <- extract_mediation(stats::glm(M ~ X, data = d, family = stats::binomial()),
                          treatment = "X", mediator = "M",
                          model_y = stats::glm(Y ~ X + M, data = d), data = d)
  testthat::local_mocked_bindings(.detect_separation = function(X, y, ...) "separated")
  r <- .pmed_model_boot(ex, 0, 1, n_boot = 5, ci_level = 0.95, seed = 1, n_sim = 500)
  expect_identical(r@diagnostics$failures[["separated"]], 5L)
  nodata <- ex
  nodata@data <- NULL
  expect_error(.pmed_model_boot(nodata, 0, 1, 5, 0.95, 1), "needs the source data")
})
