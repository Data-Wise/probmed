test_that("multilevel_designs has the documented shape and column types", {
  expect_s3_class(multilevel_designs, "data.frame")
  expect_equal(nrow(multilevel_designs), 3L)
  expect_identical(
    names(multilevel_designs),
    c(
      "study", "cohort", "clusters_sampled", "clusters_participating",
      "units_total", "mean_units_per_cluster", "cluster_response_rate",
      "restricted_use", "sampling_note", "source"
    )
  )

  expect_type(multilevel_designs$study, "character")
  expect_type(multilevel_designs$cohort, "character")
  expect_type(multilevel_designs$clusters_sampled, "integer")
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

  student <- multilevel_designs[
    multilevel_designs$study == "mediation::student",
  ]
  expect_equal(student$clusters_participating, 568L)
  expect_equal(student$units_total, 9679L)
  expect_false(student$restricted_use)
  expect_true(is.na(student$clusters_sampled))
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
