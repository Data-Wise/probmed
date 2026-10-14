# Observed-sample covariate standardization (repair tasks T9, T10).
#
# References are written from scratch (model.matrix + coef, no probmed
# helpers). With a Bernoulli outcome the half-tie P_med is linear in the arm
# means, so with a shared row r:
#   P_med = 0.5 + 0.5 * mean_r(p_t(r) - p_c(r)),
#   p_x(r) = E[ P(Y = 1 | x_value, M(x), C_r) ].
# For a Bernoulli mediator the inner expectation is a two-term sum; for a
# Gaussian mediator it is a one-dimensional integral (stats::integrate).

std_fixture <- function(n = 400, seed = 31) {
  set.seed(seed)
  X <- stats::rbinom(n, 1, 0.5)
  C <- stats::rnorm(n)
  Fc <- factor(sample(c("lo", "mid", "hi"), n, replace = TRUE))
  M <- stats::rbinom(n, 1, stats::plogis(-0.25 + 0.8 * X + 0.9 * C))
  Mg <- 0.5 * X + 0.6 * C + stats::rnorm(n)
  Y <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.3 * X + 1.2 * M + 1.5 * C +
                                           0.5 * (Fc == "hi")))
  Yg <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.3 * X + 0.9 * Mg + 1.5 * C))
  data.frame(X = X, C = C, Fc = Fc, M = M, Mg = Mg, Y = Y, Yg = Yg)
}

ref_binary_mediator <- function(fit_m, fit_y, d, x_ref = 0, x_value = 1) {
  arm <- function(x_m) {
    dm <- d; dm$X <- x_m
    pm <- stats::plogis(drop(stats::model.matrix(stats::formula(fit_m), dm) %*% stats::coef(fit_m)))
    py <- function(m) {
      dy <- d; dy$X <- x_value; dy$M <- m
      stats::plogis(drop(stats::model.matrix(stats::delete.response(stats::terms(fit_y)), dy) %*% stats::coef(fit_y)))
    }
    pm * py(1) + (1 - pm) * py(0)
  }
  0.5 + 0.5 * mean(arm(x_value) - arm(x_ref))
}

ref_gaussian_mediator <- function(fit_m, fit_y, d, x_ref = 0, x_value = 1) {
  s_m <- sqrt(sum(stats::residuals(fit_m)^2) / fit_m$df.residual)
  cm <- stats::coef(fit_m); cy <- stats::coef(fit_y)
  arm <- function(x_m) vapply(seq_len(nrow(d)), function(r) {
    mu <- cm[["(Intercept)"]] + cm[["X"]] * x_m + cm[["C"]] * d$C[r]
    lin <- cy[["(Intercept)"]] + cy[["X"]] * x_value + cy[["C"]] * d$C[r]
    stats::integrate(function(m) stats::plogis(lin + cy[["Mg"]] * m) *
                       stats::dnorm(m, mu, s_m), -Inf, Inf, rel.tol = 1e-10)$value
  }, numeric(1))
  0.5 + 0.5 * mean(arm(x_value) - arm(x_ref))
}

n_sim_big <- 4e5
mc_tol <- 4 * sqrt(0.25 / n_sim_big) # 4 Monte Carlo SEs of a half-tie mean

test_that("binary mediator and outcome with covariates match g-computation", {
  d <- std_fixture()
  fm <- stats::glm(M ~ X + C, data = d, family = stats::binomial())
  fy <- stats::glm(Y ~ X + M + C + Fc, data = d, family = stats::binomial())
  truth <- ref_binary_mediator(fm, fy, d)
  r <- pmed(Y ~ X + M + C + Fc, formula_m = M ~ X + C, data = d,
            treatment = "X", mediator = "M", family_m = stats::binomial(),
            family_y = stats::binomial(), method = "plugin",
            n_sim = n_sim_big, seed = 3)
  expect_lt(abs(r@estimate - truth), mc_tol)
})

test_that("Gaussian mediator, binary outcome with a covariate matches g-computation", {
  d <- std_fixture()
  fm <- stats::glm(Mg ~ X + C, data = d)
  fy <- stats::glm(Yg ~ X + Mg + C, data = d, family = stats::binomial())
  truth <- ref_gaussian_mediator(fm, fy, d)
  r <- pmed(Yg ~ X + Mg + C, formula_m = Mg ~ X + C, data = d,
            treatment = "X", mediator = "Mg", family_y = stats::binomial(),
            method = "plugin", n_sim = n_sim_big, seed = 4)
  expect_lt(abs(r@estimate - truth), mc_tol)
})

test_that("positive control: evaluating covariates at 0 misses the truth", {
  # The pre-repair core set every covariate to 0. That value must differ from
  # the standardized truth by well over the Monte Carlo tolerance (gap ~0.015,
  # tolerance ~0.003), or the
  # tests above could not tell the two apart.
  d <- std_fixture()
  fm <- stats::glm(M ~ X + C, data = d, family = stats::binomial())
  fy <- stats::glm(Y ~ X + M + C + Fc, data = d, family = stats::binomial())
  truth <- ref_binary_mediator(fm, fy, d)
  d0 <- d; d0$C <- 0; d0$Fc <- factor("lo", levels = levels(d$Fc))
  at_zero <- ref_binary_mediator(fm, fy, d0)
  expect_gt(abs(at_zero - truth), 3 * mc_tol)
})

test_that("Gaussian additive models are not standardized (covariates cancel)", {
  d <- std_fixture()
  d$Ygauss <- 0.3 * d$X + 0.4 * d$Mg + 0.8 * d$C + stats::rnorm(nrow(d))
  ex <- extract_mediation(stats::glm(Mg ~ X + C, data = d), treatment = "X",
                          mediator = "Mg",
                          model_y = stats::glm(Ygauss ~ X + Mg + C, data = d),
                          data = d)
  expect_false(.needs_standardization(ex))
  r <- pmed(ex, method = "plugin", n_sim = n_sim_big, seed = 5)
  closed <- stats::pnorm(ex@a_path * ex@b_path /
                           sqrt(2 * (ex@b_path^2 * ex@sigma_m^2 + ex@sigma_y^2)))
  expect_lt(abs(r@estimate - closed), mc_tol)
})

test_that("bootstraps standardize too, and the nonparametric refit handles factors", {
  d <- std_fixture(n = 300)
  args <- list(Y ~ X + M + C + Fc, formula_m = M ~ X + C, data = d,
               treatment = "X", mediator = "M", family_m = stats::binomial(),
               family_y = stats::binomial(), n_boot = 40, n_sim = 4000, seed = 6)
  for (m in c("parametric_bootstrap", "nonparametric_bootstrap")) {
    r <- do.call(pmed, c(args, method = m))
    expect_true(all(is.finite(r@boot_estimates)), label = m)
    expect_true(r@ci_lower < r@ci_upper, label = m)
  }
  # Gaussian outcome with a factor covariate: the old formula rebuild could
  # not evaluate the dummy column name "Fcmid" and failed.
  d$Yn <- 0.4 * d$Mg + 0.5 * (d$Fc == "hi") + stats::rnorm(nrow(d))
  r <- pmed(Yn ~ X + Mg + Fc, formula_m = Mg ~ X + Fc, data = d, treatment = "X",
            mediator = "Mg", method = "nonparametric_bootstrap", n_boot = 20,
            n_sim = 2000, seed = 7)
  expect_true(is.finite(r@estimate))
})

test_that("design rebuild reproduces model.matrix, including factors and transforms", {
  d <- std_fixture()
  ex <- extract_mediation(
    stats::glm(M ~ X + C + Fc, data = d, family = stats::binomial()),
    treatment = "X", mediator = "M",
    model_y = stats::glm(Y ~ X + M + log(abs(C) + 1) + Fc, data = d,
                         family = stats::binomial()),
    data = d
  )
  des <- .design_matrices(ex)
  ref_m <- stats::model.matrix(~ X + C + Fc, d)
  ref_y <- stats::model.matrix(~ X + M + log(abs(C) + 1) + Fc, d)
  expect_equal(unname(des$m), unname(ref_m[, , drop = FALSE]), ignore_attr = TRUE)
  expect_equal(unname(des$y), unname(ref_y[, , drop = FALSE]), ignore_attr = TRUE)
  expect_identical(colnames(des$y), colnames(ref_y))
})
