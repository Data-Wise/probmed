# Composite-null diagnostics for method = "mbco" (Phase 8, T3).
# At a = b = 0 the profile LR at p* = 0.5 is the minimum of the two branch
# statistics, so the chi-square(1) interval is conservative (about 0.9975).

mbco_null_data <- function(a, b, n = 400, seed = 1) {
  set.seed(seed)
  d <- data.frame(X = stats::rbinom(n, 1, 0.5), C = stats::rnorm(n))
  d$M <- a * d$X + 0.3 * d$C + stats::rnorm(n)
  d$Y <- 0.3 * d$X + b * d$M + 0.2 * d$C + stats::rnorm(n)
  d
}
fit_mbco <- function(d) {
  pmed(Y ~ X + M + C, formula_m = M ~ X + C, data = d, treatment = "X",
       mediator = "M", method = "mbco")
}

test_that("branch LR statistics equal the nested-model likelihood ratios (known answer)", {
  d <- mbco_null_data(0.3, 0.2)
  ex <- medfit::extract_mediation(
    stats::lm(M ~ X + C, d), model_y = stats::lm(Y ~ X + M + C, d),
    treatment = "X", mediator = "M", data = d)
  prep <- .mbco_prep(ex, 0, 1)
  lr <- .mbco_branch_lr(prep)
  ref_a <- 2 * (as.numeric(stats::logLik(stats::lm(M ~ X + C, d))) -
                as.numeric(stats::logLik(stats::lm(M ~ C, d))))
  ref_b <- 2 * (as.numeric(stats::logLik(stats::lm(Y ~ X + M + C, d))) -
                as.numeric(stats::logLik(stats::lm(Y ~ X + C, d))))
  expect_equal(unname(lr[["a0"]]), ref_a, tolerance = 1e-8)
  expect_equal(unname(lr[["b0"]]), ref_b, tolerance = 1e-8)
})

test_that("near_double_null is TRUE at a = b = 0 and FALSE when either path is strong", {
  r0 <- fit_mbco(mbco_null_data(0, 0, seed = 3))
  expect_true(isTRUE(r0@diagnostics$pmed$near_double_null))
  expect_true(all(r0@diagnostics$pmed$branch_lr < stats::qchisq(0.95, 1)))

  r_a <- fit_mbco(mbco_null_data(0.8, 0, seed = 3)) # single null b = 0, a strong
  expect_false(r_a@diagnostics$pmed$near_double_null)
  r_b <- fit_mbco(mbco_null_data(0, 0.6, seed = 3)) # single null a = 0, b strong
  expect_false(r_b@diagnostics$pmed$near_double_null)
  r_ab <- fit_mbco(mbco_null_data(0.8, 0.6, seed = 3))
  expect_false(r_ab@diagnostics$pmed$near_double_null)
})

test_that("the flag uses the requested level's cutoff", {
  d <- mbco_null_data(0.12, 0.12, n = 400, seed = 11)
  r <- fit_mbco(d)
  lr <- r@diagnostics$pmed$branch_lr
  expect_identical(r@diagnostics$pmed$near_double_null,
                   all(lr < stats::qchisq(0.95, 1)))
})

test_that("print notes conservativeness only when near the double null", {
  r0 <- fit_mbco(mbco_null_data(0, 0, seed = 3))
  r1 <- fit_mbco(mbco_null_data(0.8, 0.6, seed = 3))
  expect_match(paste(utils::capture.output(print(r0)), collapse = "\n"),
               "conservative", fixed = TRUE)
  expect_no_match(paste(utils::capture.output(print(r1)), collapse = "\n"),
                  "conservative", fixed = TRUE)
})

test_that(".mbco_lr_at is ~0 at the estimate, matches the interval cutoff at its endpoints, and is 0.5-composite", {
  d <- mbco_null_data(0.4, 0.3, seed = 5)
  r <- fit_mbco(d)
  prep <- .mbco_prep(r@source_extract, 0, 1)
  expect_lt(abs(.mbco_lr_at(prep, r@estimate)), 1e-6)
  crit <- stats::qchisq(0.95, 1)
  expect_equal(.mbco_lr_at(prep, r@ci_lower), crit, tolerance = 1e-4)
  expect_equal(.mbco_lr_at(prep, r@ci_upper), crit, tolerance = 1e-4)
  # at p = 0.5 the statistic is the minimum of the two branch statistics
  expect_equal(.mbco_lr_at(prep, 0.5), min(.mbco_branch_lr(prep)), tolerance = 1e-8)
})
