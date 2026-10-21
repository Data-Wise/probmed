## Coverage for the core PmedResult print/summary/plot methods (methods-print.R)
## and the PmedResult validator branches (classes.R). These are the package's
## primary user-facing output; every other estimator's print method is tested,
## but the flagship pmed() result's was not (methods-print.R was at 0%).

.pp_data <- function(n = 200, seed = 1) {
  set.seed(seed)
  data <- data.frame(X = stats::rnorm(n), C = stats::rnorm(n))
  data$M <- 0.5 * data$X + 0.3 * data$C + stats::rnorm(n)
  data$Y <- 0.4 * data$M + 0.2 * data$X + 0.2 * data$C + stats::rnorm(n)
  data
}

.pp_fit <- function(method, n_boot = 200L, seed = 1) {
  pmed(Y ~ X + M + C, formula_m = M ~ X + C, data = .pp_data(seed = seed),
       treatment = "X", mediator = "M", method = method, n_boot = n_boot)
}

test_that("print() covers the bootstrap result (CI + IE + n_boot branches)", {
  r <- .pp_fit("parametric_bootstrap")
  out <- capture.output(print(r))
  expect_true(any(grepl("P_med", out)))
  expect_true(any(grepl("CI: \\[", out)))          # ci_lower not NA -> CI branch
  expect_true(any(grepl("Indirect Effect", out)))  # ie_estimate not NA branch
  expect_true(any(grepl("Bootstrap samples", out)))# n_boot not NA branch
  expect_true(any(grepl("Interpretation", out)))
})

test_that("print() covers the plugin result (no-CI branch)", {
  r <- .pp_fit("plugin")
  out <- capture.output(print(r))
  expect_true(any(grepl("P_med", out)))
  # plugin has no CI: the "% CI: [" line must be absent
  expect_false(any(grepl("CI: \\[", out)))
})

test_that("summary() covers both the bootstrap-distribution and plain branches", {
  out_boot <- capture.output(summary(.pp_fit("parametric_bootstrap")))
  expect_true(any(grepl("Bootstrap Distribution", out_boot)))
  expect_true(any(grepl("Standard Error", out_boot)))

  out_plug <- capture.output(summary(.pp_fit("plugin")))
  # plugin: no boot_estimates -> distribution block absent, but still prints source
  expect_false(any(grepl("Bootstrap Distribution", out_plug)))
  expect_true(any(grepl("Sample size", out_plug)))
})

test_that("plot() draws for a bootstrap result and errors without a distribution", {
  r_boot <- .pp_fit("parametric_bootstrap")
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  expect_invisible(plot(r_boot))   # returns invisible(x)
  grDevices::dev.off()
  unlink(tf)

  # plugin has no bootstrap distribution -> the guard must stop()
  expect_error(plot(.pp_fit("plugin")), "No bootstrap distribution")
})

test_that("PmedResult validator flags out-of-range estimate and bad ci bounds", {
  # estimate outside [0,1] -> warning branch (classes.R)
  expect_warning(
    PmedResult(estimate = 1.5, ci_lower = NA_real_, ci_upper = NA_real_,
               ci_level = 0.95, method = "plugin", x_ref = 0, x_value = 1,
               source_extract = NULL, converged = TRUE),
    "outside \\[0,1\\]"
  )
  # ci_lower > ci_upper -> validator error branch
  expect_error(
    PmedResult(estimate = 0.5, ci_lower = 0.8, ci_upper = 0.2,
               ci_level = 0.95, method = "plugin", x_ref = 0, x_value = 1,
               source_extract = NULL, converged = TRUE),
    "ci_lower must be"
  )
  # ci_level outside (0,1) -> validator error branch
  expect_error(
    PmedResult(estimate = 0.5, ci_lower = NA_real_, ci_upper = NA_real_,
               ci_level = 1.5, method = "plugin", x_ref = 0, x_value = 1,
               source_extract = NULL, converged = TRUE),
    "ci_level must be"
  )
})

test_that("null-coalescing operator returns fallback only on NULL", {
  `%||%` <- probmed:::`%||%`
  expect_equal(NULL %||% 3, 3)
  expect_equal(5 %||% 3, 5)
  expect_equal(NA %||% 3, NA)  # NA is not NULL -> returns NA
})
