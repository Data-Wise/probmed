test_that("incr_sensitivity runs end-to-end on an incr_pmed fit", {
  set.seed(1)
  n <- 800
  C <- rnorm(n)
  A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
  M <- 0.6 * A + 0.4 * C + rnorm(n)
  Y <- 0.5 * A + 0.7 * M + 0.3 * C + rnorm(n)
  fit <- incr_pmed(data.frame(A, M, Y, C), deltas = c(0.5, 1, 2))
  sens <- incr_sensitivity(fit, threshold = 0)
  expect_s3_class(sens, "data.frame")
  expect_equal(nrow(sens), 3L)
  expect_equal(sens$tipping, -fit@curve$med)
  expect_true(all(is.finite(sens$tipping)))
})

test_that("incr_sensitivity returns one row per delta with closed-form tipping bias", {
  curve <- data.frame(
    delta = c(0.5, 1, 2),
    dir = c(0.20, 0.25, 0.30),
    med = c(0.10, 0.15, 0.05),
    tot = c(0.30, 0.40, 0.35),
    Pmed = c(0.10, 0.15, 0.05) / c(0.30, 0.40, 0.35),
    se = c(0.02, 0.02, 0.02), lo = NA, hi = NA
  )
  res <- probmed:::.incr_sensitivity_from_curve(curve, threshold = 0)
  expect_equal(nrow(res), 3L)
  expect_equal(res$delta, curve$delta)
  expect_equal(res$tipping, -curve$med)                 # b = -med
  expect_equal(res$tipping_threshold, 0 * curve$tot - curve$med)
  expect_equal(res$Pmed, curve$med / curve$tot)
})

test_that("incr_sensitivity tipping_threshold is distinct from tipping when threshold != 0", {
  curve <- data.frame(
    delta = c(0.5, 1, 2),
    dir = c(0.20, 0.25, 0.30),
    med = c(0.10, 0.15, 0.05),
    tot = c(0.30, 0.40, 0.35),
    Pmed = c(0.10, 0.15, 0.05) / c(0.30, 0.40, 0.35),
    se = c(0.02, 0.02, 0.02), lo = NA, hi = NA
  )
  res <- probmed:::.incr_sensitivity_from_curve(curve, threshold = 0.5)
  expect_equal(res$tipping, -curve$med)
  expect_equal(res$tipping_threshold, 0.5 * curve$tot - curve$med)
})
