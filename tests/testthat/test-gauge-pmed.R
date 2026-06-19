## Tests for ward_residual() / GaugePmedResult (gauge-calibrated P_med)
.gp_gen <- function(n, tint, binaryY, seed = 1) {
  set.seed(seed)
  C <- rnorm(n); A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
  M <- 0.6 * A + 0.4 * C + rnorm(n)
  lin <- 0.5 * A + 0.7 * M + tint * A * M + 0.3 * C
  Y <- if (binaryY) rbinom(n, 1, plogis(lin)) else lin + rnorm(n)
  data.frame(A, M, Y, C)
}

test_that("ward_residual returns a GaugePmedResult with finite estimates", {
  r <- ward_residual(.gp_gen(1500, 0.0, FALSE))
  expect_true(S7::S7_inherits(r, GaugePmedResult))
  expect_true(is.finite(r@p_med)); expect_true(is.finite(r@W))
})

test_that("identity R = OE - IDE - IIE holds numerically", {
  r <- ward_residual(.gp_gen(2000, 0.5, FALSE))
  expect_equal(r@R, r@OE - r@IDE - r@IIE, tolerance = 1e-8)
})

test_that("no interaction => W ~ 0 and CI covers 0", {
  r <- ward_residual(.gp_gen(3000, 0.0, FALSE))
  expect_lt(abs(r@W), 0.1)
  expect_lte(r@W_ci[1], 0); expect_gte(r@W_ci[2], 0)
})

test_that("A*M interaction => W > 0 with CI excluding 0", {
  r <- ward_residual(.gp_gen(4000, 0.9, TRUE))
  expect_gt(r@W, 0.1)
  expect_gt(r@W_ci[1], 0)
  expect_lt(r@W_p, 0.01)
})

test_that("Fieller set is computed and typed", {
  r <- ward_residual(.gp_gen(3000, 0.5, FALSE), fieller = TRUE)
  expect_true(r@fieller_type %in% c("bounded", "exclusive-unbounded", "all-real", "empty"))
  expect_length(r@p_med_fieller, 2)
})

test_that("fieller=FALSE leaves the Fieller fields empty", {
  r <- ward_residual(.gp_gen(1500, 0.3, FALSE), fieller = FALSE)
  expect_length(r@p_med_fieller, 0)
  expect_true(is.na(r@fieller_type))
})

test_that("multiple covariates are supported", {
  d <- .gp_gen(1500, 0.3, FALSE); d$C2 <- rnorm(nrow(d))
  r <- ward_residual(d, covars = c("C", "C2"))
  expect_true(is.finite(r@W))
})
