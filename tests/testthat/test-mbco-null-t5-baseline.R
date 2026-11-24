t5_baseline_dir <- testthat::test_path("..", "..", "inst", "sim")
if (!dir.exists(t5_baseline_dir)) {
  t5_baseline_dir <- system.file("sim", package = "probmed", mustWork = TRUE)
}

test_that("T5 reconstructs branch-specific t magnitudes and selects by t order", {
  env <- new.env()
  sys.source(file.path(t5_baseline_dir, "calibration", "calib_mbco_null_t5_baseline.R"), env)
  # First row is frozen task 1, replication 558. Its LR and t order differ.
  fixture <- data.frame(
    N = c(100L, 1000L), cov = c(FALSE, TRUE),
    lr_a0 = c(0.3649889, 0.25), lr_b0 = c(0.3664547, 1)
  )
  branches <- env$reconstruct_t5_branches(fixture)
  expect_equal(fixture$N * log1p(branches$t_a^2 / c(98, 997)), fixture$lr_a0)
  expect_equal(fixture$N * log1p(branches$t_b^2 / c(97, 996)), fixture$lr_b0)
  expect_identical(branches$a_first, c(FALSE, TRUE))
  expect_identical(branches$lr_disagrees, c(TRUE, FALSE))
  expect_equal(branches$branch_df, c(97, 997))
  expect_identical(env$t5_disposition(c("INCONCLUSIVE", "FAIL")), "FAIL")
  expect_identical(env$t5_disposition(c("PASS", "INCONCLUSIVE")), "INCONCLUSIVE")
  expect_identical(env$t5_disposition(c("PASS", "PASS")), "PASS")
})

test_that("T5 reproduces the frozen T4 baselines and branch-order discrepancy", {
  lattice <- file.path(t5_baseline_dir, "results", "calibration", "mbco_null")
  skip_if_not(file.exists(file.path(lattice, "replications.csv.gz")))
  env <- new.env()
  sys.source(file.path(t5_baseline_dir, "calibration", "calib_mbco_null_score.R"), env)
  sys.source(file.path(t5_baseline_dir, "calibration", "calib_mbco_null_t5_baseline.R"), env)
  reps <- utils::read.csv(file.path(lattice, "replications.csv.gz"))
  saved <- utils::read.csv(file.path(lattice, "scores.csv"))
  audit <- env$audit_t5_baseline(reps, saved, env$score_rules)

  expect_equal(nrow(audit$scores), 400)
  expect_true(all(audit$scores$reps == 2000))
  expect_equal(audit$branch_disagreements, 157)
  expect_true(all(audit$gates$overall == "FAIL"))
  chi2_no_cov <- subset(audit$scores, rule == "chi2" & N == 1000 & !cov)
  expect_equal(min(chi2_no_cov$coverage), 0.9440)
  expect_equal(sum(!is.na(audit$scores$origin_gate)), 16)
  saved$ci_lower[1] <- saved$ci_lower[1] - 1e-6
  expect_error(env$audit_t5_baseline(reps, saved, env$score_rules), "ci_lower")
})

test_that("MC21 t uses the direct checkerboard test only on the null", {
  env <- new.env()
  sys.source(file.path(t5_baseline_dir, "calibration", "calib_mbco_null_t5_baseline.R"), env)
  sys.source(file.path(t5_baseline_dir, "calibration", "calib_mbco_null_t5_score.R"), env)
  # Folded t probabilities 0.01 and 0.08 are in different bands.
  lr <- function(u, df) 100 * log1p(stats::qt((1 + u) / 2, df)^2 / df)
  fixture <- data.frame(
    N = 100L, cov = FALSE, h_a = c(0, 0, 1, 1), h_b = 1,
    lr_a0 = lr(0.01, 98), lr_b0 = lr(c(0.01, 0.08, 0.01, 0.08), 97),
    lr_true = c(lr(0.01, 98), lr(0.01, 98), 1, 4)
  )
  accept <- env$t5_prototype_acceptance(fixture)

  expect_identical(accept[["MC21_t_05"]], c(FALSE, TRUE, TRUE, FALSE))
  expect_false("MC21_t_all" %in% names(accept))
})

test_that("T5 prototype scores regenerate from the frozen lattice", {
  lattice <- file.path(t5_baseline_dir, "results", "calibration", "mbco_null")
  skip_if_not(file.exists(file.path(lattice, "replications.csv.gz")))
  env <- new.env()
  for (script in c("calib_mbco_null_score.R", "calib_mbco_null_t5_baseline.R",
                   "calib_mbco_null_t5_score.R")) {
    sys.source(file.path(t5_baseline_dir, "calibration", script), env)
  }
  reps <- utils::read.csv(file.path(lattice, "replications.csv.gz"))
  baseline <- utils::read.csv(file.path(lattice, "scores.csv"))
  saved <- utils::read.csv(file.path(lattice, "t5_prototype_scores.csv"))
  audit <- env$audit_t5_prototypes(reps, baseline, env$score_rules)

  expect_equal(audit$scores, saved, tolerance = 1e-12)
  expect_equal(nrow(saved), 1500)
  expect_true(all(table(saved$rule) == 100L))
  expect_true(all(saved$reps == 2000L))
  expect_false("MC21_t_all" %in% saved$rule)
  expect_false("PASS" %in% audit$gates$overall)
  mc21 <- subset(saved, rule == "MC21_t_05" & h_a == 2 & h_b == 0 & N == 1000 & !cov)
  expect_equal(mc21$accepted, 1876)
  expect_identical(mc21$minimum_gate, "INCONCLUSIVE")
})
