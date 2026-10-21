# Known-answer checks for the composite-null scoring script (Phase 8, T4).
score_path <- testthat::test_path("..", "..", "inst", "sim", "calibration", "calib_mbco_null_score.R")

test_that("scoring reproduces min-of-two chi-square theory and flags an undercovering rule", {
  skip_if_not(file.exists(score_path))
  env <- new.env()
  sys.source(score_path, envir = env)
  set.seed(1)
  n <- 4e4
  dn <- data.frame(h_a = 0, h_b = 0, N = 1000, cov = FALSE,
                   lr_a0 = stats::rchisq(n, 1), lr_b0 = stats::rchisq(n, 1))
  dn$lr_true <- pmin(dn$lr_a0, dn$lr_b0)            # double null: min of branches
  sn <- data.frame(h_a = 4, h_b = 0, N = 1000, cov = FALSE,
                   lr_a0 = stats::rchisq(n, 1, ncp = 16), lr_b0 = stats::rchisq(n, 1))
  sn$lr_true <- sn$lr_b0                            # single null b = 0
  sc <- env$score_rules(rbind(dn, sn))
  get <- function(rule, h) sc$coverage[sc$rule == rule & sc$h_a == h]
  expect_equal(get("chi2", 0), 1 - 0.05^2, tolerance = 0.002)    # 0.9975 theory
  expect_equal(get("chi2", 4), 0.95, tolerance = 0.005)
  expect_lt(get("switch_logN", 0), get("chi2", 0))                # less conservative...
  expect_lt(get("switch_logN", 4), 0.945)                         # ...but undercovers here
  s <- env$summarize_rules(sc)
  expect_false(s$pass[s$rule == "chi2"])                          # fails the origin criterion
})
