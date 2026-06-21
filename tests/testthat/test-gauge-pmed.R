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

## ---- P1 Fix-path: reps (repeated cross-fitting) + bootstrap se ----

test_that("se_method and reps are recorded on the result", {
  r <- ward_residual(.gp_gen(800, 0, FALSE), se_method = "bootstrap", reps = 3L, B = 80L)
  expect_equal(r@se_method, "bootstrap")
  expect_equal(r@reps, 3L)
})

test_that("se_method defaults to analytic with reps = 1", {
  r <- ward_residual(.gp_gen(1500, 0, FALSE))
  expect_equal(r@se_method, "analytic")
  expect_equal(r@reps, 1L)
})

test_that("se_method='bootstrap' yields a different, positive W_se than analytic", {
  d <- .gp_gen(800, 0, FALSE)
  ra <- ward_residual(d, se_method = "analytic")
  rb <- ward_residual(d, se_method = "bootstrap", B = 150L)
  expect_gt(rb@W_se, 0)
  expect_false(isTRUE(all.equal(ra@W_se, rb@W_se)))
})

test_that("bootstrap W_ci is the percentile interval, not the symmetric Wald form", {
  d <- .gp_gen(800, 0, FALSE)
  rb <- ward_residual(d, se_method = "bootstrap", B = 200L)
  ## percentile interval contains the point and is NOT W +/- z*se (analytic-style)
  expect_lte(rb@W_ci[1], rb@W); expect_gte(rb@W_ci[2], rb@W)
  wald <- c(rb@W - qnorm(0.975) * rb@W_se, rb@W + qnorm(0.975) * rb@W_se)
  expect_false(isTRUE(all.equal(rb@W_ci, wald, tolerance = 1e-6)))
})

test_that("bootstrap intervals are reproducible under a fixed seed", {
  d <- .gp_gen(800, 0, FALSE)
  r1 <- ward_residual(d, se_method = "bootstrap", B = 120L, seed = 7L)
  r2 <- ward_residual(d, se_method = "bootstrap", B = 120L, seed = 7L)
  expect_equal(r1@W_ci, r2@W_ci); expect_equal(r1@p_med_ci, r2@p_med_ci)
})

test_that("analytic W_ci remains the symmetric Wald interval", {
  rb <- ward_residual(.gp_gen(1500, 0.3, FALSE))   # default analytic
  expect_equal(rb@W_ci, c(rb@W - qnorm(0.975) * rb@W_se,
                          rb@W + qnorm(0.975) * rb@W_se), tolerance = 1e-8)
})

test_that("reps>1 averages corner influence over independent fold draws", {
  d <- .gp_gen(2000, 0.5, FALSE)
  r1 <- ward_residual(d, reps = 1L)
  r4 <- ward_residual(d, reps = 4L)
  expect_equal(r4@reps, 4L)
  expect_lt(abs(r4@W - r1@W), 0.05)          # reps-averaged point stays close
  expect_true(is.finite(r4@W_se) && r4@W_se > 0)
})

test_that("bootstrap se preserves the identity R = OE - IDE - IIE", {
  r <- ward_residual(.gp_gen(1500, 0.5, FALSE), se_method = "bootstrap", B = 100L)
  expect_equal(r@R, r@OE - r@IDE - r@IIE, tolerance = 1e-8)
})

test_that("print labels the interval by construction (Wald vs percentile)", {
  d <- .gp_gen(800, 0.5, FALSE)
  ra <- ward_residual(d, se_method = "analytic")
  rb <- ward_residual(d, se_method = "bootstrap", B = 80L)
  oa <- paste(capture.output(print(ra)), collapse = "\n")
  ob <- paste(capture.output(print(rb)), collapse = "\n")
  expect_match(oa, "Wald")
  expect_no_match(oa, "percentile")
  expect_match(ob, "percentile")
  expect_no_match(ob, "Wald")
})
