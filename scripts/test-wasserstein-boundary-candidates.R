#!/usr/bin/env Rscript
devtools::load_all(quiet = TRUE)
source("scripts/wasserstein-boundary-candidates.R")
testthat::test_that("damped correction preserves the target geometry", {
  v <- matrix(c(0.0001, -0.00002, -0.00002, 0.0002), 2)
  d <- c(0.12, 0.10)
  p <- wb_ratio(d)
  testthat::expect_equal(wb_correct(rev(d), v[2:1, 2:1]), 1 - wb_correct(d, v))
  testthat::expect_equal(wb_correct(3 * d, 9 * v), wb_correct(d, v))
  testthat::expect_equal(wb_correct(d, v, k = 1e12), p, tolerance = 1e-10)
  for (x in c(0, 1e-12, 0.001, 0.1, -0.1)) {
    a <- wb_correct(c(x, 0.1), v)
    testthat::expect_true(is.finite(a) && a >= 0 && a <= 1)
  }
  testthat::expect_equal(wb_correct(c(0, 0.1), v), 0)
  testthat::expect_equal(wb_correct(c(0.1, 0), v), 1)
  testthat::expect_error(wb_correct(c(0, 0), v), "undefined")
  tiny <- v * 1e-5
  hessian <- numDeriv::hessian(wb_ratio, d)
  expected <- p - sum(hessian * tiny) / 2
  testthat::expect_equal(wb_correct(d, tiny), expected, tolerance = 1e-10)
  # Independent finite-difference Hessian must agree with the derived expression.
  a <- c(0.5 / d[1], -0.5 / d[2])
  analytic <- p * (1 - p) * ((1 - 2 * p) * tcrossprod(a) +
    diag(c(-0.5 / d[1]^2, 0.5 / d[2]^2)))
  testthat::expect_equal(hessian, analytic, tolerance = 1e-7)
})
testthat::test_that("independent binary prototype matches production inference", {
  for (seed in c(19L, 37L, 20281186L)) {
    set.seed(seed)
    C <- stats::rnorm(500)
    A <- stats::rbinom(500, 1, stats::plogis(0.3 * C))
    M <- 0.7 * A + 0.5 * C + stats::rnorm(500)
    Y <- stats::rbinom(500, 1, stats::plogis(-0.5 + 0.6 * A + M + 0.3 * C))
    d <- data.frame(A, M, Y, C)
    prototype <- wb_binary(d)
    production <- probmed:::.wasserstein_observed(d, binY = TRUE)
    testthat::expect_equal(prototype$raw, production$pmedW, tolerance = 1e-9)
    testthat::expect_equal(prototype$se, production$se, tolerance = 1e-8)
  }
})
