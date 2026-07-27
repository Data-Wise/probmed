# data-raw/multilevel_designs.R
#
# Build `multilevel_designs`: publicly documented design metadata for a small
# set of clustered (multilevel) study designs, so users planning a clustered
# mediation analysis can see realistic numbers of clusters.
#
# Every value below is transcribed from the public source named in the row's
# `source` field. Nothing here is estimated, modeled, or derived from any
# analysis, with the single exception of `mean_units_per_cluster`, which is
# computed as units_total / clusters_participating (see below).
#
# Re-create the shipped object with:
#   source("data-raw/multilevel_designs.R")

multilevel_designs <- data.frame(
  study = c(
    "ECLS-K",
    "ECLS-K",
    "mediation::student"
  ),
  cohort = c(
    "1998-99",
    "2010-11",
    NA_character_
  ),
  # Canonical "sampled" count. For ECLS-K:2011 this is the original school
  # sample drawn before substitution; the eligible and substituted counts are
  # recorded in `sampling_note` rather than in separate columns.
  clusters_sampled = c(
    1280L,
    1352L,
    NA_integer_
  ),
  clusters_participating = c(
    940L,
    970L,
    568L
  ),
  units_total = c(
    21260L,
    18000L,
    9679L
  ),
  cluster_response_rate = c(
    0.74,
    NA_real_,
    NA_real_
  ),
  restricted_use = c(
    TRUE,
    TRUE,
    FALSE
  ),
  sampling_note = c(
    paste(
      "Weighted school response rate as reported by NCES.",
      "Child count is the base-year kindergarten sample."
    ),
    paste(
      "Original school sample before substitution was 1352;",
      "1264 of those were eligible for the fall collection;",
      "93 substitute schools were added, giving 1446 total recruitment",
      "attempts. Participating-school and child counts are reported as",
      "approximate in the user's manual. No school response rate is given",
      "in the cited source, so the rate is NA."
    ),
    paste(
      "Public teaching data set shipped with the mediation R package.",
      "No sampled-school count or response rate is documented, so both",
      "are NA."
    )
  ),
  source = c(
    "NCES Handbook of Survey Methods, https://nces.ed.gov/statprog/handbook/ecls_k_dataquality.asp",
    "User's Manual for the ECLS-K:2011 Kindergarten Data File (NCES 2015-074), https://nces.ed.gov/pubs2015/2015074.pdf",
    "mediation R package (`data(student, package = \"mediation\")`)"
  ),
  stringsAsFactors = FALSE
)

# Derived, not sourced: the denominator is the *participating* cluster count,
# not the sampled count.
multilevel_designs$mean_units_per_cluster <-
  multilevel_designs$units_total / multilevel_designs$clusters_participating

# Fix column order (derived column inserted after units_total).
multilevel_designs <- multilevel_designs[, c(
  "study",
  "cohort",
  "clusters_sampled",
  "clusters_participating",
  "units_total",
  "mean_units_per_cluster",
  "cluster_response_rate",
  "restricted_use",
  "sampling_note",
  "source"
)]

rownames(multilevel_designs) <- NULL

# Guard against non-ASCII sneaking into a shipped data set.
chr_cols <- vapply(multilevel_designs, is.character, logical(1))
stopifnot(
  !any(grepl("[^\x01-\x7f]", unlist(multilevel_designs[chr_cols]), useBytes = TRUE))
)

usethis::use_data(multilevel_designs, overwrite = TRUE)
