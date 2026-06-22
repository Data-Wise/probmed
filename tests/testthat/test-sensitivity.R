test_that("no-confounding case reduces to the point estimate", {
  s <- pmed_sensitivity(indirect = 0.30, total = 0.50)
  ## b = 0 is on the default symmetric grid (midpoint).
  i0 <- which(abs(s@bias_grid) < 1e-12)
  expect_length(i0, 1L)
  expect_equal(s@p_med_bias[i0], s@p_med)
  expect_equal(s@p_med, 0.30 / 0.50)
})

test_that("known-input numeric checks for bias propagation", {
  s <- pmed_sensitivity(indirect = 0.30, total = 0.50, threshold = 0)
  ## P_med(b) = (0.3 + b) / 0.5
  expect_equal(s@p_med_bias, (0.30 + s@bias_grid) / 0.50)
  ## tipping bias zeroing the indirect effect: b = -indirect.
  expect_equal(s@tipping_indirect, -0.30)
  ## tipping bias for P_med -> 0 : b = 0*total - indirect = -indirect.
  expect_equal(s@tipping_threshold, -0.30)
})

test_that("threshold tipping bias solves (indirect + b)/total = threshold", {
  s <- pmed_sensitivity(indirect = 0.20, total = 0.80, threshold = 0.10)
  b <- s@tipping_threshold
  expect_equal((0.20 + b) / 0.80, 0.10)
  ## explicit value: b = 0.10*0.80 - 0.20 = -0.12
  expect_equal(b, -0.12)
})

test_that("custom bias grid is honored and propagated", {
  g <- c(-0.1, 0, 0.1)
  s <- pmed_sensitivity(indirect = 0.40, total = 1.00, bias_grid = g)
  expect_equal(s@bias_grid, g)
  expect_equal(s@p_med_bias, c(0.30, 0.40, 0.50))
})

test_that("extracts components from a GaugePmedResult", {
  fit <- GaugePmedResult(
    p_med = 0.6, p_med_ci = c(0.4, 0.8),
    W = 0.05, W_ci = c(-0.1, 0.2), W_se = 0.07, W_p = 0.5,
    OE = 0.5, IDE = 0.2, IIE = 0.3, R = 0,
    theta = c(0.5, 0.2, 0.3, 0), method = "test", n = 100L, ci_level = 0.95
  )
  s <- pmed_sensitivity(fit, threshold = 0)
  expect_equal(s@indirect, 0.3)
  expect_equal(s@total, 0.5)
  expect_equal(s@p_med, 0.6)
  expect_equal(s@tipping_indirect, -0.3)
})

test_that("direct indirect/total override an object", {
  fit <- GaugePmedResult(
    p_med = 0.6, p_med_ci = c(0.4, 0.8),
    W = 0, W_ci = c(0, 0), W_se = 0, W_p = 1,
    OE = 0.5, IDE = 0.2, IIE = 0.3, R = 0,
    theta = c(0.5, 0.2, 0.3, 0), method = "test", n = 100L, ci_level = 0.95
  )
  s <- pmed_sensitivity(fit, indirect = 0.10, total = 0.40)
  expect_equal(s@indirect, 0.10)
  expect_equal(s@total, 0.40)
})

test_that("zero or non-finite denominator errors", {
  expect_error(pmed_sensitivity(indirect = 0.3, total = 0), "undefined")
  expect_error(pmed_sensitivity(indirect = 0.3, total = NA_real_), "undefined")
})

test_that("missing components raise an informative error", {
  expect_error(pmed_sensitivity(object = list(foo = 1)),
               "indirect/total")
})

test_that("print method runs without error", {
  s <- pmed_sensitivity(indirect = 0.3, total = 0.5, threshold = 0,
                        scale = "Delta")
  expect_output(print(s), "sensitivity")
  expect_output(print(s), "Tipping")
})

test_that("E-value summary is guarded by EValue availability", {
  if (requireNamespace("EValue", quietly = TRUE)) {
    s <- pmed_sensitivity(indirect = 0.3, total = 0.5,
                          evalue = TRUE, se = 0.1)
    expect_true(length(s@evalue) >= 1L)
  } else {
    expect_warning(
      pmed_sensitivity(indirect = 0.3, total = 0.5, evalue = TRUE),
      "EValue"
    )
  }
})
