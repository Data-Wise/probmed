# Common random numbers in .pmed_core_simple() (task T21).

g <- stats::gaussian()
bl <- stats::binomial()
exact_gauss <- function(a, b, sm, sy) stats::pnorm(a * b / sqrt(2 * (b^2 * sm^2 + sy^2)))

test_that("the same u gives the same P_med; no u leaves the RNG path unchanged", {
  set.seed(1)
  u <- probmed:::.crn_uniforms(2000)
  f <- function(uu) probmed:::.pmed_core_simple(
    a = 0.5, b = 0.4, c_prime = 0.3, x_ref = 0, x_value = 1, n_sim = 2000,
    family_y = g, family_m = g, u = uu)
  expect_identical(f(u), f(u))
  set.seed(7); r1 <- f(NULL)
  set.seed(7)
  r2 <- probmed:::.pmed_core_simple(a = 0.5, b = 0.4, c_prime = 0.3, x_ref = 0,
                                    x_value = 1, n_sim = 2000, family_y = g, family_m = g)
  expect_identical(r1, r2)
})

test_that("CRN Gaussian branch hits the closed form (known answer)", {
  set.seed(2)
  u <- probmed:::.crn_uniforms(4e5)
  est <- probmed:::.pmed_core_simple(
    a = 0.5, b = 0.4, c_prime = 0, x_ref = 0, x_value = 1, sigma_m = 1.2,
    sigma_y = 0.8, n_sim = 4e5, family_y = g, family_m = g, u = u)
  expect_equal(est, exact_gauss(0.5, 0.4, 1.2, 0.8), tolerance = 3e-3)
})

test_that("CRN Bernoulli branches use u < p and match the RNG path in expectation", {
  set.seed(3)
  args <- list(a = 0.8, b = 1.1, c_prime = 0.2, x_ref = 0, x_value = 1, n_sim = 4e5,
               family_y = bl, family_m = bl, i_m = -0.3, i_y = -0.5)
  crn <- do.call(probmed:::.pmed_core_simple, c(args, list(u = probmed:::.crn_uniforms(4e5))))
  rng <- do.call(probmed:::.pmed_core_simple, args)
  expect_equal(crn, rng, tolerance = 4e-3)
  expect_gt(crn, 0.5)
})

test_that("CRN standardization rows stay inside the analysis rows", {
  u <- cbind(c(0, 0.5, 1 - 1e-12, 1), matrix(0.5, 4, 4))
  est <- probmed:::.pmed_core_simple(
    a = 0.5, b = 0.4, c_prime = 0, x_ref = 0, x_value = 1, n_sim = 4,
    family_y = bl, family_m = bl, lp_m0 = c(-1, 1), lp_y0 = c(0, 0), u = u)
  expect_true(is.finite(est))
})

test_that("CRN removes Monte Carlo spread between bootstrap-style evaluations", {
  # Positive control: identical coefficients, many evaluations. Fresh draws
  # vary by about 0.005 at n_sim = 1e4; one shared u gives zero spread.
  set.seed(4)
  u <- probmed:::.crn_uniforms(1e4)
  f <- function(uu) probmed:::.pmed_core_simple(
    a = 0.5, b = 0.4, c_prime = 0, x_ref = 0, x_value = 1, n_sim = 1e4,
    family_y = g, family_m = g, u = uu)
  fresh <- replicate(50, f(NULL))
  shared <- replicate(50, f(u))
  expect_gt(stats::sd(fresh), 0.002)
  expect_equal(stats::sd(shared), 0)
})

test_that("model_bootstrap interval narrows to the exact-evaluation width", {
  set.seed(5)
  n <- 1000
  d <- data.frame(X = stats::rbinom(n, 1, 0.5), C = stats::rbinom(n, 1, 0.5))
  d$M <- 0.5 * d$X + 0.4 * d$C + stats::rnorm(n)
  d$Y <- 0.3 * d$X + 0.4 * d$M + 0.4 * d$C + stats::rnorm(n)
  r <- suppressWarnings(pmed(Y ~ X + M + C, formula_m = M ~ X + C, data = d,
                             treatment = "X", mediator = "M", method = "model_bootstrap",
                             n_boot = 400, n_sim = 1e4, seed = 11))
  # T20 (exact P_med per draw) median width at this design is 0.0296; the
  # fresh-draw interval was 0.0354.
  w <- r@ci_upper - r@ci_lower
  expect_lt(w, 0.034)
  expect_gt(w, 0.022)
  fm <- stats::lm(M ~ X + C, d); fy <- stats::lm(Y ~ X + M + C, d)
  ex <- exact_gauss(stats::coef(fm)[["X"]], stats::coef(fy)[["M"]],
                    stats::sigma(fm), stats::sigma(fy))
  # Point estimate uses 10 * n_sim pairs: MC sd about 0.0016.
  expect_equal(r@estimate, ex, tolerance = 0.006)
})

test_that("draws are centered on the point estimate (control variate)", {
  # Positive control for the T22 pilot failure (coverage 0.79 with an
  # uncentered shared u). Deliberately biased uniforms make P_u far from
  # P_med; the anchor P_u(fit) must absorb that bias so the interval still
  # brackets the point estimate. Binary M and Y with a covariate, so the
  # standardization and Bernoulli branches both run.
  set.seed(6)
  n <- 1000
  d <- data.frame(X = stats::rbinom(n, 1, 0.5), C = stats::rbinom(n, 1, 0.5))
  d$M <- stats::rbinom(n, 1, stats::plogis(-0.25 + 0.5 * d$X + 0.4 * d$C))
  d$Y <- stats::rbinom(n, 1, stats::plogis(-0.5 + 0.3 * d$X + 0.4 * d$M + 0.4 * d$C))
  local_mocked_bindings(.crn_uniforms = function(n_sim) {
    u <- matrix(stats::runif(5L * n_sim), ncol = 5L)
    u[, 2] <- u[, 2]^3 # M(x) almost always 1: P_u biased upward
    u
  })
  r <- suppressWarnings(pmed(Y ~ X + M + C, formula_m = M ~ X + C, data = d,
                             treatment = "X", mediator = "M",
                             family_m = stats::binomial(), family_y = stats::binomial(),
                             method = "model_bootstrap", n_boot = 200, n_sim = 1e4,
                             seed = 3))
  expect_true(r@ci_lower < r@estimate && r@estimate < r@ci_upper)
  expect_lt(r@ci_upper - r@ci_lower, 0.05)
})
