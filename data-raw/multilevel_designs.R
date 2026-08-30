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
  # "probability sample" = a designed sample with a documented frame;
  # "teaching data" = a convenience data set shipped for instruction.
  design_type = c(
    "probability sample",
    "probability sample",
    "teaching data"
  ),
  # Cluster counts along the recruitment chain. The four columns have
  # DIFFERENT bases, which is why they are kept separate rather than folded
  # into a single "sampled" count:
  #   clusters_sampled       original sample drawn (before substitution)
  #   clusters_eligible      sampled clusters found eligible for collection
  #   clusters_recruited     total recruitment attempts, including substitutes
  #   clusters_participating clusters that actually took part (may include
  #                          substitutes, so it is not comparable to
  #                          clusters_sampled)
  # A cell is NA_integer_ whenever the cited source documents no such count.
  # Do not fill it in from another source or by inference.
  clusters_sampled = c(
    1280L,
    1352L,
    NA_integer_
  ),
  clusters_eligible = c(
    NA_integer_,
    1264L,
    NA_integer_
  ),
  clusters_recruited = c(
    NA_integer_,
    1446L,
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
      "attempts. The participating count includes substitutes, so",
      "participating / sampled is not a response rate.",
      "Participating-school and child counts are reported as approximate",
      "in the user's manual. No school response rate is given in the",
      "cited source, so the rate is NA."
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
  "design_type",
  "clusters_sampled",
  "clusters_eligible",
  "clusters_recruited",
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

# version = 2 keeps the serialization readable by R >= 2.10; the v3 default
# would force a user-visible `Depends: R (>= 3.5)` for a three-row table.
usethis::use_data(multilevel_designs, overwrite = TRUE, version = 2)
