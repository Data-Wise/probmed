## testthat tests for gauge_pmed.R (run standalone or move to probmed/tests/testthat/)
library(testthat); library(S7)
source("gauge_pmed.R")
expit <- function(x) 1/(1+exp(-x))
gen <- function(n, tint, binaryY, seed=1) {
  set.seed(seed)
  C <- rnorm(n); A <- rbinom(n, 1, expit(-0.2 + 0.8*C))
  M <- 0.6*A + 0.4*C + rnorm(n)
  lin <- 0.5*A + 0.7*M + tint*A*M + 0.3*C
  Y <- if (binaryY) rbinom(n, 1, expit(lin)) else lin + rnorm(n)
  data.frame(A, M, Y, C)
}

test_that("returns GaugePmedResult with finite estimates", {
  r <- ward_residual(gen(1500, 0.0, FALSE))
  expect_s3_class(r, "probmed::GaugePmedResult")  # S7 objects also carry S3 class tag
  expect_true(is.finite(r@p_med)); expect_true(is.finite(r@W))
})

test_that("no interaction => W ~ 0 (CI covers 0)", {
  r <- ward_residual(gen(3000, 0.0, FALSE))
  expect_lt(abs(r@W), 0.1)
  expect_lte(r@W_ci[1], 0); expect_gte(r@W_ci[2], 0)
})

test_that("A*M interaction => W significantly > 0", {
  r <- ward_residual(gen(4000, 0.9, TRUE))
  expect_gt(r@W, 0.1)
  expect_gt(r@W_ci[1], 0)            # CI excludes 0
  expect_lt(r@W_p, 0.01)
})

test_that("identity R = OE - IDE - IIE holds numerically", {
  r <- ward_residual(gen(2000, 0.5, FALSE))
  expect_equal(r@R, r@OE - r@IDE - r@IIE, tolerance = 1e-8)
})
