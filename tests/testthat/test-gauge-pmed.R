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

## ---- non-binary-coded (two-level) exposure: a0 / a1 contrast ----

test_that("a0/a1 default reproduces the 0/1 binary result", {
  d <- .gp_gen(1200, 0.5, FALSE)
  r_default <- ward_residual(d, seed = 3L)
  r_explicit <- ward_residual(d, a0 = 0, a1 = 1, seed = 3L)
  expect_equal(r_default@W, r_explicit@W)
  expect_equal(r_default@p_med, r_explicit@p_med)
})

test_that("factor / non-0-1 coding matches the recoded 0/1 fit", {
  d <- .gp_gen(1200, 0.5, FALSE)
  d_fac <- d; d_fac$A <- factor(ifelse(d$A == 1, "trt", "ctrl"))
  d_num <- d; d_num$A <- ifelse(d$A == 1, 2L, 1L)
  r01 <- ward_residual(d, seed = 5L)
  rfac <- ward_residual(d_fac, a0 = "ctrl", a1 = "trt", seed = 5L)
  rnum <- ward_residual(d_num, a0 = 1L, a1 = 2L, seed = 5L)
  expect_equal(rfac@W, r01@W, tolerance = 1e-8)
  expect_equal(rnum@p_med, r01@p_med, tolerance = 1e-8)
})

test_that("reversing a0/a1 flips the sign of W but preserves |W|", {
  d <- .gp_gen(1200, 0.5, FALSE)
  r_fwd <- ward_residual(d, a0 = 0, a1 = 1, seed = 9L)
  r_rev <- ward_residual(d, a0 = 1, a1 = 0, seed = 9L)
  ## W = R/OE: the mixed 2nd difference R is invariant to the swap, OE flips sign,
  ## so W flips sign while |W| (the non-decomposability magnitude) is preserved.
  expect_equal(r_fwd@W, -r_rev@W, tolerance = 1e-6)
  expect_equal(abs(r_fwd@W), abs(r_rev@W), tolerance = 1e-6)
})

test_that("A with >2 levels errors (no silent restricted-estimand)", {
  d <- .gp_gen(1200, 0.5, FALSE)
  d$A[1:50] <- 2L                       # introduce a third level
  expect_error(ward_residual(d, a0 = 0, a1 = 1), "level")
})

test_that("a0/a1 not present in A errors", {
  d <- .gp_gen(800, 0.5, FALSE)
  expect_error(ward_residual(d, a0 = 0, a1 = 5), "not found|level")
})

# ---- Issue #11: weak-identification flag (A2) + regularity guard (A1) --------

# near-null DGM: A affects neither M nor Y, so OE -> 0 (non-regular ratio).
.gp_null <- function(n, seed = 1) {
  set.seed(seed)
  C <- rnorm(n); A <- rbinom(n, 1, 0.5)
  M <- 0.4 * C + rnorm(n); Y <- 0.7 * M + 0.3 * C + rnorm(n)
  data.frame(A, M, Y, C)
}

test_that("A1: well-identified OE gives oe_regular = TRUE, high oe_snr", {
  r <- ward_residual(.gp_gen(2000, 0.9, FALSE))
  expect_true(is.finite(r@oe_snr) && r@oe_snr > 2)
  expect_true(isTRUE(r@oe_regular))
})

test_that("A1: near-singular OE flags oe_regular = FALSE and warns under bootstrap", {
  ws <- character()
  r <- withCallingHandlers(
    ward_residual(.gp_null(800), se_method = "bootstrap", B = 60L),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  expect_true(any(grepl("near-singular|non-regular", ws, ignore.case = TRUE)))
  expect_false(isTRUE(r@oe_regular))
  expect_lt(r@oe_snr, 2)
})

test_that("A1: near-singular OE does NOT warn under analytic (no bootstrap CI reported)", {
  expect_no_warning(ward_residual(.gp_null(800), se_method = "analytic", fieller = FALSE))
})

test_that("A2: weak_id is NA under analytic (no percentile interval to compare)", {
  r <- ward_residual(.gp_gen(1500, 0.5, FALSE), se_method = "analytic")
  expect_true(is.na(r@weak_id))
  expect_true(is.na(r@weak_id_ratio))
})

test_that("A2: strongly-identified case does NOT trip the weak-ID flag", {
  # well-identified: percentile is ~2x wider than Wald (anti-conservatism), below 3x
  r <- ward_residual(.gp_gen(3000, 0.9, FALSE), se_method = "bootstrap", B = 200L)
  expect_true(is.logical(r@weak_id) && !is.na(r@weak_id))
  expect_true(is.finite(r@weak_id_ratio) && r@weak_id_ratio > 0)
  expect_false(r@weak_id)
  expect_lt(r@weak_id_ratio, 3)
})

test_that("A2: near-singular case trips the weak-ID flag (percentile >= 3x Wald)", {
  suppressWarnings(
    r <- ward_residual(.gp_null(800), se_method = "bootstrap", B = 200L)
  )
  expect_true(r@weak_id_ratio >= 3)
  expect_true(isTRUE(r@weak_id))
})

test_that("gate fields are present on the result", {
  r <- ward_residual(.gp_gen(800, 0.5, FALSE))
  expect_length(r@W_ci_wald, 2)
  for (f in c("weak_id", "weak_id_ratio", "oe_snr", "oe_regular"))
    expect_true(f %in% S7::prop_names(r))
})

test_that("gate thresholds are tunable via arguments", {
  d <- .gp_gen(1200, 0.9, FALSE)
  r_default <- suppressWarnings(ward_residual(d, se_method = "bootstrap", B = 100L))
  expect_false(r_default@weak_id)  # default threshold=3, well-ID case
  r_strict <- suppressWarnings(
    ward_residual(d, se_method = "bootstrap", B = 100L, weak_id_ratio_threshold = 1)
  )
  expect_true(r_strict@weak_id)  # threshold=1 must trip on any width inflation
  r_lenient <- ward_residual(d, oe_snr_threshold = 0)
  expect_true(isTRUE(r_lenient@oe_regular))  # threshold=0: any finite oe_snr passes
  # non-tautological direction: a threshold ABOVE a well-identified case's oe_snr
  # must flip oe_regular to FALSE, proving the comparison actually uses the arg.
  r_overstrict <- ward_residual(d, oe_snr_threshold = 100)
  expect_true(r_overstrict@oe_snr < 100)
  expect_false(isTRUE(r_overstrict@oe_regular))
})

test_that("print() surfaces the weak-ID and near-singular-OE flags with the ACTUAL threshold used", {
  r_weak <- suppressWarnings(
    ward_residual(.gp_gen(1200, 0.9, FALSE), se_method = "bootstrap", B = 100L,
                  weak_id_ratio_threshold = 1)
  )
  expect_output(print(r_weak), "weak-ID.*wider than Wald")
  # regression guard: print() must report the threshold actually applied (1),
  # not a hardcoded default (3) -- the bug this test was written to catch.
  expect_output(print(r_weak), ">= 1x")
  expect_no_match(capture.output(print(r_weak)), ">= 3x", all = FALSE)

  r_nonreg <- suppressWarnings(
    ward_residual(.gp_null(800), se_method = "bootstrap", B = 60L, oe_snr_threshold = 5)
  )
  expect_output(print(r_nonreg), "near-singular OE.*non-regular")
  expect_output(print(r_nonreg), "< 5")
})

test_that("weak_id and oe_regular stay NA (not FALSE) in degenerate zero-width/zero-se cases", {
  # constructing GaugePmedResult directly to exercise the NA-preserving branches
  # without needing to engineer a genuinely zero-width bootstrap CI or zero se(OE)
  # from real data (astronomically unlikely, but the code path must not crash or
  # silently report a confident FALSE for an undefined comparison).
  r <- GaugePmedResult(
    p_med = 0.5, p_med_ci = c(0.4, 0.6), W = 0.1, W_ci = c(0.05, 0.15),
    W_se = 0.02, W_p = 0.1, W_ci_wald = c(0.08, 0.08),  # zero-width Wald
    weak_id = if (is.na(NA_real_)) NA else FALSE, weak_id_ratio = NA_real_,
    oe_snr = NA_real_, oe_regular = NA,
    OE = 0.5, IDE = 0.1, IIE = 0.3, R = 0.1, theta = c(0, 0, 0, 0),
    method = "onestep-crossfit", n = 100L, ci_level = 0.95, se_method = "bootstrap"
  )
  expect_true(is.na(r@weak_id))
  expect_true(is.na(r@oe_regular))
  # print() must not error when these flags are NA (isTRUE/isFALSE guard, not bare if)
  expect_no_error(print(r))
})

## ---- corner-EIF weight fix (2026-08-22) and the reps > 1 se construction ----

test_that("corner EIF weights are per-row: corner means are invariant to row order", {
  ## Regression test for the ifelse() weight bug (shipped 10b41aa, fixed 2026-08-22):
  ## pa(z) <- ifelse(z == 1, p1, 1 - p1) returned p1[1] -- the FIRST row's propensity
  ## -- for every row, so the corner means depended on which row came first. The
  ## full-sample fit consumes no RNG, so a row permutation must leave them unchanged.
  d <- .gp_gen(600, 0.5, FALSE)
  m1 <- colMeans(.corner_fit_full(d, FALSE, "C")$phi)
  m2 <- colMeans(.corner_fit_full(d[rev(seq_len(nrow(d))), ], FALSE, "C")$phi)
  expect_equal(m1, m2, tolerance = 1e-10)
  expect_gt(stats::sd(.corner_fit_full(d, FALSE, "C")$g), 0.05)  # weights do vary by row
})

test_that("reps = 1 analytic path: values pinned after the weight fix (regression guard)", {
  set.seed(1); n <- 800; C <- rnorm(n)
  A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C)); M <- 0.6 * A + 0.4 * C + rnorm(n)
  Y <- 0.5 * A + 0.7 * M + 0.8 * A * M + 0.3 * C + rnorm(n)
  r1 <- ward_residual(data.frame(A, M, Y, C))
  ## pre-fix (dev @ b3ea449): W 0.3753, W_se 0.0947, p_med CI [0.210, 0.332]
  expect_equal(r1@W, 0.4028370648, tolerance = 1e-8)
  expect_equal(r1@W_se, 0.0503393240, tolerance = 1e-8)
  expect_equal(r1@p_med_ci, c(0.1922297676, 0.3376918320), tolerance = 1e-8)
})

test_that("reps > 1 point estimate pinned after the weight fix; se finite and positive", {
  set.seed(1); n <- 800; C <- rnorm(n)
  A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C)); M <- 0.6 * A + 0.4 * C + rnorm(n)
  Y <- 0.5 * A + 0.7 * M + 0.8 * A * M + 0.3 * C + rnorm(n)
  r4 <- ward_residual(data.frame(A, M, Y, C), reps = 4L)
  expect_equal(r4@W, 0.3973313916, tolerance = 1e-8)   # pre-fix: 0.4146
  expect_true(is.finite(r4@W_se) && r4@W_se > 0)
})

test_that(".corner_fit_full matches the cross-fit corner means at large n and consumes no RNG", {
  d <- .gp_gen(4000, 0.5, FALSE)
  set.seed(3); cf <- colMeans(.corner_fit(d, 5L, FALSE, "C")$phi)
  s0 <- .Random.seed
  full <- .corner_fit_full(d, FALSE, "C")
  expect_identical(.Random.seed, s0)
  expect_equal(dim(full$phi), c(4000L, 4L))
  expect_equal(colnames(full$phi), c("11", "10", "01", "00"))
  expect_equal(unname(colMeans(full$phi)), unname(cf), tolerance = 0.05)
})

test_that("reps > 1 analytic se = averaged cross-fit IF se (+) residual fold Monte-Carlo term", {
  d <- .gp_gen(800, 0.5, FALSE)
  r <- ward_residual(d, reps = 4L, seed = 11L)
  ## rebuild by hand: the same per-rep fold draws (seed + r), the averaged phi, its
  ## IF with the averaged W, plus var(W_reps)/reps.
  d2 <- d; d2$A <- as.integer(d2$A == 1)
  W_reps <- numeric(4); phis <- vector("list", 4)
  for (k in 1:4) {
    set.seed(11L + k); phis[[k]] <- .corner_fit(d2, 5L, FALSE, "C")$phi
    t <- colMeans(phis[[k]]); oe <- t["11"] - t["00"]
    W_reps[k] <- unname((oe - (t["10"] - t["00"]) - (t["01"] - t["00"])) / oe)
  }
  pbar <- Reduce(`+`, phis) / 4
  pOE <- pbar[, "11"] - pbar[, "00"]
  pR  <- pOE - (pbar[, "10"] - pbar[, "00"]) - (pbar[, "01"] - pbar[, "00"])
  se_cf <- stats::sd((pR - r@W * pOE) / mean(pOE)) / sqrt(800)
  expect_equal(r@W_se, sqrt(se_cf^2 + stats::var(W_reps) / 4), tolerance = 1e-10)
})

test_that("reps > 1 analytic se works for a binary outcome", {
  r <- ward_residual(.gp_gen(800, 0.5, TRUE), reps = 2L)
  expect_true(is.finite(r@W_se) && r@W_se > 0)
  expect_true(is.finite(r@p_med_ci[1]) && is.finite(r@p_med_ci[2]))
})
