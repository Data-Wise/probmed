test_that("multilevel_designs has the documented shape and column types", {
  expect_s3_class(multilevel_designs, "data.frame")
  expect_equal(nrow(multilevel_designs), 3L)
  expect_identical(
    names(multilevel_designs),
    c(
      "study", "cohort", "design_type", "clusters_sampled",
      "clusters_eligible", "clusters_recruited", "clusters_participating",
      "units_total", "mean_units_per_cluster", "cluster_response_rate",
      "restricted_use", "sampling_note", "source"
    )
  )

  expect_type(multilevel_designs$study, "character")
  expect_type(multilevel_designs$cohort, "character")
  expect_type(multilevel_designs$design_type, "character")
  expect_type(multilevel_designs$clusters_sampled, "integer")
  expect_type(multilevel_designs$clusters_eligible, "integer")
  expect_type(multilevel_designs$clusters_recruited, "integer")
  expect_type(multilevel_designs$clusters_participating, "integer")
  expect_type(multilevel_designs$units_total, "integer")
  expect_type(multilevel_designs$mean_units_per_cluster, "double")
  expect_type(multilevel_designs$cluster_response_rate, "double")
  expect_type(multilevel_designs$restricted_use, "logical")
  expect_type(multilevel_designs$sampling_note, "character")
  expect_type(multilevel_designs$source, "character")
})

test_that("multilevel_designs carries the documented source values", {
  eclsk98 <- multilevel_designs[
    multilevel_designs$study == "ECLS-K" &
      multilevel_designs$cohort == "1998-99",
  ]
  expect_equal(eclsk98$clusters_sampled, 1280L)
  expect_equal(eclsk98$clusters_participating, 940L)
  expect_equal(eclsk98$units_total, 21260L)
  expect_equal(eclsk98$cluster_response_rate, 0.74)
  expect_true(eclsk98$restricted_use)
  expect_equal(eclsk98$design_type, "probability sample")
  # The 1998-99 source documents no eligible or recruited count.
  expect_true(is.na(eclsk98$clusters_eligible))
  expect_true(is.na(eclsk98$clusters_recruited))

  student <- multilevel_designs[
    multilevel_designs$study == "mediation::student",
  ]
  expect_equal(student$clusters_participating, 568L)
  expect_equal(student$units_total, 9679L)
  expect_false(student$restricted_use)
  expect_true(is.na(student$clusters_sampled))
  expect_true(is.na(student$clusters_eligible))
  expect_true(is.na(student$clusters_recruited))
  expect_true(is.na(student$cluster_response_rate))
  expect_equal(student$design_type, "teaching data")
})

test_that("the ECLS-K:2011 row carries the user's-manual recruitment chain", {
  # NCES 2015-074 (base-year user's manual). An NCES summary page implies
  # roughly 1,310 participating schools; the manual gives roughly 970. This
  # table carries the manual's figure. A "correction" to 1310 must fail here.
  eclsk11 <- multilevel_designs[
    multilevel_designs$study == "ECLS-K" &
      multilevel_designs$cohort == "2010-11",
  ]
  expect_equal(nrow(eclsk11), 1L)
  expect_equal(eclsk11$design_type, "probability sample")
  expect_equal(eclsk11$clusters_sampled, 1352L)
  expect_equal(eclsk11$clusters_eligible, 1264L)
  expect_equal(eclsk11$clusters_recruited, 1446L)
  expect_equal(eclsk11$clusters_participating, 970L)
  expect_equal(eclsk11$units_total, 18000L)
  expect_true(eclsk11$restricted_use)

  # The chain is ordered: sampled >= eligible, recruited >= eligible. No
  # arithmetic relation among the stages is asserted beyond what the source
  # states, since the source reports the substitute count and the total
  # attempts separately.
  expect_gte(eclsk11$clusters_sampled, eclsk11$clusters_eligible)
  expect_gte(eclsk11$clusters_recruited, eclsk11$clusters_eligible)

  # No response rate is reported by the source, and participating / sampled
  # must not be offered as one (different bases).
  expect_true(is.na(eclsk11$cluster_response_rate))
  expect_match(eclsk11$sampling_note, "not a response rate", fixed = TRUE)
})

test_that("cluster counts are NA only where undocumented, never invented", {
  # Exactly one row (ECLS-K:2011) documents eligible and recruited counts.
  expect_equal(sum(!is.na(multilevel_designs$clusters_eligible)), 1L)
  expect_equal(sum(!is.na(multilevel_designs$clusters_recruited)), 1L)
  expect_equal(sum(!is.na(multilevel_designs$clusters_sampled)), 2L)
  expect_false(anyNA(multilevel_designs$clusters_participating))
})

test_that("mean_units_per_cluster is units_total / clusters_participating", {
  expect_equal(
    multilevel_designs$mean_units_per_cluster,
    multilevel_designs$units_total / multilevel_designs$clusters_participating
  )
})

test_that("multilevel_designs is ASCII-only", {
  chr_cols <- vapply(multilevel_designs, is.character, logical(1))
  values <- unlist(multilevel_designs[chr_cols], use.names = FALSE)
  expect_false(any(grepl("[^\x01-\x7f]", values, useBytes = TRUE)))
})
