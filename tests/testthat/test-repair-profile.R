repair_prep <- function() {
  d <- generate_mediation_data(n = 250, a = 0.5, b = 0.4, seed = 142)
  ex <- extract_mediation(stats::lm(M ~ X, d), treatment = "X", mediator = "M",
                          model_y = stats::lm(Y ~ X + M, d), data = d)
  .mbco_prep(ex, 0, 1)
}

test_that("cached residual quadratics match direct least squares", {
  prep <- repair_prep()
  for (coefficient in c(-2, 0, 0.5, 2)) {
    expect_equal(.mbco_sse(prep$quad_m, coefficient),
      sum(stats::residuals(stats::lm.fit(prep$Dm, prep$M - coefficient * prep$X))^2),
      tolerance = 1e-10)
    expect_equal(.mbco_sse(prep$quad_y, coefficient),
      sum(stats::residuals(stats::lm.fit(prep$Dy, prep$Y - coefficient * prep$M))^2),
      tolerance = 1e-10)
  }
})

test_that("signed profiling reconstructs the target in both basins", {
  prep <- repair_prep()
  for (delta in c(1, -2)) {
    prep$delta <- delta
    for (q in c(-0.2, 0.2)) {
      fit <- .mbco_fit_pmed(prep, q)
      expect_true(fit$converged)
      expect_equal(length(fit$basins), 2)
      for (basin in fit$basins) {
        p <- basin$parameters
        actual <- p$a * delta * p$b / sqrt(2 * (p$b^2 * p$Vm + p$Vy))
        expect_equal(actual, q, tolerance = 1e-8)
        expect_gt(p$Vm, 0)
        expect_gt(p$Vy, 0)
      }
      expect_lte(fit$ll, prep$ll_free + 1e-6)
    }
  }
})

test_that("signed profiles agree with an independent full Gaussian likelihood", {
  skip_if_not_installed("nloptr")
  d <- generate_mediation_data(n = 180, a = 0.5, b = 0.4, seed = 155)
  set.seed(156)
  d$C <- stats::rnorm(nrow(d))
  d$M <- d$M + 0.4 * d$C
  d$Y <- d$Y + 0.3 * d$C
  ex <- extract_mediation(stats::lm(M ~ X + C, d), treatment = "X", mediator = "M",
                          model_y = stats::lm(Y ~ X + M + C, d), data = d)
  for (delta in c(1, -2)) for (target in c(0.4, 0.6)) {
    prep <- .mbco_prep(ex, 0, delta)
    actual <- .mbco_ll_constrained_pmed(prep, stats::qnorm(target), 1)
    reference <- repair_full_oracle(d, target, delta)
    refined <- repair_full_oracle(d, target, delta, scales = c(0.25, 4))
    expect_true(is.finite(actual))
    expect_lt(abs(refined - reference) / nrow(d), 1e-6)
    expect_lt(abs(actual - reference) / nrow(d), 1e-6)
  }
})

test_that("optimizer exhaustion is never accepted as a constrained fit", {
  prep <- repair_prep()
  # nloptr status 5 = maxeval reached (exhaustion), never a success code.
  testthat::local_mocked_bindings(
    nloptr = function(x0, ...) list(status = 5L, objective = 10, solution = x0),
    .package = "nloptr"
  )
  expect_true(is.na(.mbco_ll_constrained_pmed(prep, 0.2, 1)))
  expect_true(is.na(.mbco_ll_constrained_ie(prep, 0.1)))
})

test_that("indirect-effect profiles agree with the independent full likelihood", {
  skip_if_not_installed("nloptr")
  d <- generate_mediation_data(n = 180, a = 0.5, b = 0.4, seed = 155)
  set.seed(156)
  d$C <- stats::rnorm(nrow(d))
  ex <- extract_mediation(stats::lm(M ~ X + C, d), treatment = "X", mediator = "M",
                          model_y = stats::lm(Y ~ X + M + C, d), data = d)
  prep <- .mbco_prep(ex, 0, 1)
  for (target in c(-0.1, 0.1)) {
    actual <- .mbco_ll_constrained_ie(prep, target)
    reference <- repair_full_oracle(d, target, estimand = "ie")
    expect_true(is.finite(actual))
    expect_lt(abs(actual - reference) / nrow(d), 1e-6)
  }
})

test_that("ML point and verified endpoints agree with independent references", {
  skip_if_not_installed("nloptr")
  d <- generate_mediation_data(n = 180, a = 0.5, b = 0.4, seed = 155)
  set.seed(156)
  d$C <- stats::rnorm(nrow(d))
  d$M <- d$M + 0.4 * d$C
  d$Y <- d$Y + 0.3 * d$C
  fm <- stats::lm(M ~ X + C, d); fy <- stats::lm(Y ~ X + M + C, d)
  ex <- extract_mediation(fm, treatment = "X", mediator = "M", model_y = fy, data = d)
  free <- as.numeric(stats::logLik(fm)) + as.numeric(stats::logLik(fy))
  for (delta in c(1, -2)) {
    point <- as.numeric(stats::pnorm(stats::coef(fm)["X"] * delta * stats::coef(fy)["M"] /
      sqrt(2 * (stats::coef(fy)["M"]^2 * mean(stats::residuals(fm)^2) +
                  mean(stats::residuals(fy)^2)))))
    result <- pmed(ex, method = "mbco", x_ref = 0, x_value = delta)
    expect_equal(result@estimate, point, tolerance = 1e-12)
    expect_true(result@converged)
    excess <- function(p) 2 * (free - repair_full_oracle(d, p, delta)) - stats::qchisq(.95, 1)
    reference <- vapply(c(-1, 1), function(direction) {
      inner <- point
      for (k in 1:20) {
        outer <- point + direction * .01 * k
        if (excess(outer) > 0) return(stats::uniroot(excess, sort(c(inner, outer)), tol = 1e-7)$root)
        inner <- outer
      }
      stop("Independent endpoint bracket not found")
    }, numeric(1))
    expect_lt(abs(result@ci_lower - reference[1]), 1e-4)
    expect_lt(abs(result@ci_upper - reference[2]), 1e-4)
  }
})

test_that("verified inversion detects roots, truncation, failure and re-entry", {
  root <- .mbco_invert_verified(0.5, function(x) (x - 0.5)^2 - 0.04,
                                c(0, 1), 0.01)
  expect_identical(root$status, "ok")
  expect_equal(root$bounds, c(lower = 0.3, upper = 0.7), tolerance = 1e-5)
  expect_identical(root$endpoint_status, c(lower = "root", upper = "root"))
  flat <- .mbco_invert_verified(0.5, function(x) -1, c(0, 1), 0.01)
  expect_identical(flat$status, "truncated")
  expect_false(flat$verified)
  failed <- .mbco_invert_verified(0.5, function(x) if (x > 0.8) NA_real_ else -1,
                                  c(0, 1), 0.01)
  expect_identical(failed$status, "failed")
  expect_true(all(is.na(failed$bounds)))
  # Two disjoint accepted intervals: [.2,.4] and [.6,.8].
  split <- .mbco_invert_verified(0.3,
    function(x) (x - 0.2) * (x - 0.4) * (x - 0.6) * (x - 0.8), c(0, 1), 0.01)
  expect_identical(split$status, "disconnected")
  expect_true(all(is.na(split$bounds)))
  expect_equal(nrow(split$components), 2)
})

test_that("MBCO exposes honest failure diagnostics", {
  d <- generate_mediation_data(n = 100, seed = 151)
  ex <- extract_mediation(stats::lm(M ~ X, d), treatment = "X", mediator = "M",
                          model_y = stats::lm(Y ~ X + M, d), data = d)
  free <- .mbco_prep(ex, 0, 1)$ll_free
  testthat::local_mocked_bindings(.mbco_ll_constrained_pmed = function(prep, qstar, ...) {
    if (abs(qstar) > 0.08) NA_real_ else free
  })
  result <- pmed(ex, method = "mbco")
  expect_false(result@converged)
  expect_true(is.na(result@ci_lower))
  expect_identical(result@diagnostics$pmed$status, "failed")
  expect_output(print(result), "Interval status.*failed")
})

test_that("near-null P_med profiles reach the b ~ 0 branch of each sign basin", {
  skip_if_not_installed("nloptr")
  # Regression (inst/sim/mbco_na_rate.R, cell strong/n=100/covariate, seed 18):
  # without a start aimed at the b ~ 0 branch, the optimizer stayed on the a ~ 0
  # branch for p just above 0.5, the profile jumped by ~6.7 between adjacent
  # points, root verification failed, and the interval came back NA.
  set.seed(1200018)
  n <- 100
  X <- stats::rnorm(n)
  C1 <- stats::rnorm(n)
  M <- 0.5 * X + 0.3 * C1 + stats::rnorm(n)
  Y <- 0.2 * X + 0.4 * M + 0.3 * C1 + stats::rnorm(n)
  d <- data.frame(X = X, M = M, Y = Y, C = C1)
  ex <- extract_mediation(stats::lm(M ~ X + C, d), treatment = "X", mediator = "M",
                          model_y = stats::lm(Y ~ X + M + C, d), data = d)
  prep <- .mbco_prep(ex, 0, 1)
  for (target in c(0.5003, 0.5006)) {
    actual <- .mbco_ll_constrained_pmed(prep, stats::qnorm(target), 1)
    reference <- repair_full_oracle(d, target, 1)
    expect_lt(abs(actual - reference) / n, 1e-6)
  }
  res <- pmed(Y ~ X + M + C, formula_m = M ~ X + C, data = d,
              treatment = "X", mediator = "M", method = "mbco")
  expect_identical(res@diagnostics$pmed$status, "ok")
  crit <- stats::qchisq(0.95, 1)
  lower_ll <- repair_full_oracle(d, res@ci_lower, 1)
  expect_lt(abs(-2 * (lower_ll - prep$ll_free) - crit), 1e-3)
})

test_that("verified inversion reports each failure path explicitly", {
  # Invalid inputs: estimate outside the domain, empty domain.
  expect_identical(.mbco_invert_verified(2, function(x) x - 1, c(0, 1), 0.1)$reason,
                   "invalid_domain_or_point")
  expect_identical(.mbco_invert_verified(0.5, function(x) x, c(1, 0), 0.1)$reason,
                   "invalid_domain_or_point")
  # The estimate itself is rejected.
  expect_identical(.mbco_invert_verified(0.5, function(x) 1, c(0, 1), 0.1)$reason,
                   "estimate_rejected")
  # A jump discontinuity has no root: verification must refuse the bracket.
  jump <- function(x) if (x < 0.3 || x > 0.7) 1 else -1
  res <- .mbco_invert_verified(0.5, jump, c(0, 1), 0.1)
  expect_identical(res$reason, "root_verification")
  expect_false(res$verified)
  expect_true(all(is.na(res$bounds)))
  # Positive control: a continuous version of the same shape verifies.
  smooth <- function(x) abs(x - 0.5) - 0.2
  ok <- .mbco_invert_verified(0.5, smooth, c(0, 1), 0.1)
  expect_identical(ok$status, "ok")
  expect_equal(unname(ok$bounds), c(0.3, 0.7), tolerance = 1e-8)
})

test_that(".mbco_optim falls back to Subplex when BOBYQA does not succeed", {
  calls <- character(0)
  testthat::local_mocked_bindings(
    nloptr = function(x0, eval_f, opts, ...) {
      calls <<- c(calls, opts$algorithm)
      if (opts$algorithm == "NLOPT_LN_BOBYQA") {
        list(status = 5L, objective = 10, solution = x0)       # maxeval: not a success
      } else {
        list(status = 4L, objective = 1.5, solution = x0 + 1)  # xtol reached: success
      }
    },
    .package = "nloptr"
  )
  o <- .mbco_optim(c(0, 0), function(p) sum(p^2))
  expect_identical(calls, c("NLOPT_LN_BOBYQA", "NLOPT_LN_SBPLX"))
  expect_identical(o$convergence, 0L)
  expect_equal(o$par, c(1, 1))
  expect_equal(o$value, 1.5)
})

test_that(".mbco_nlopt turns an nloptr error into NULL, and a real problem converges", {
  # a real nloptr error (the objective throws), so the path is hit without mocking
  expect_null(.mbco_nlopt(c(0, 0), function(p) stop("boom"), "NLOPT_LN_BOBYQA", maxeval = 10L))
})

test_that(".mbco_optim minimizes a known quadratic (real nloptr)", {
  o <- .mbco_optim(c(3, -2), function(p) sum((p - c(1, 2))^2) + 7)
  expect_identical(o$convergence, 0L)
  expect_equal(o$par, c(1, 2), tolerance = 1e-5)
  expect_equal(o$value, 7, tolerance = 1e-8)
})
