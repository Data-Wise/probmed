test_that("explicit plugin seed is honored regardless of caller RNG state", {
  set.seed(1)
  d <- data.frame(x = stats::rbinom(200, 1, 0.5))
  d$m <- 0.5 * d$x + stats::rnorm(200)
  d$y <- 0.4 * d$m + stats::rnorm(200)
  run <- function() pmed(y ~ x + m, formula_m = m ~ x, data = d,
                        treatment = "x", mediator = "m", method = "plugin",
                        seed = 42, n_sim = 1000)
  set.seed(101)
  first <- run()@estimate
  set.seed(202)
  expect_identical(run()@estimate, first)
})

test_that("both exact-null likelihoods maximize over both reduced models", {
  for (branch in c("a", "b")) {
    set.seed(17)
    x <- rep(c(0, 1), 100)
    c <- stats::rnorm(200)
    em <- stats::residuals(stats::lm(stats::rnorm(200) ~ x + c))
    m <- if (branch == "a") em else 1.5 * x + em
    ey <- stats::residuals(stats::lm(stats::rnorm(200) ~ x + c + m))
    d <- data.frame(x = x, c = c, m = m,
                    y = 0.3 * x + 0.4 * c + if (branch == "a") 1.5 * m + ey else ey)
    ex <- extract_mediation(stats::lm(m ~ x + c, d), treatment = "x",
                            mediator = "m", model_y = stats::lm(y ~ x + m + c, d),
                            data = d)
    prep <- .mbco_prep(ex, 0, 1)
    a0 <- as.numeric(stats::logLik(stats::lm(m ~ c, d))) +
      as.numeric(stats::logLik(stats::lm(y ~ x + m + c, d)))
    b0 <- as.numeric(stats::logLik(stats::lm(m ~ x + c, d))) +
      as.numeric(stats::logLik(stats::lm(y ~ x + c, d)))
    expect_gt(if (branch == "a") a0 - b0 else b0 - a0, 1)
    expected <- max(a0, b0)
    expect_equal(.mbco_ll_constrained_ie(prep, 0), expected, tolerance = 1e-8)
    expect_equal(.mbco_ll_constrained_pmed(prep, 0, 1), expected, tolerance = 1e-8)
  }
})
