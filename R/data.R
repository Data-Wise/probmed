#' Design metadata for candidate multilevel studies
#'
#' Publicly documented design metadata for a small set of clustered
#' (multilevel) studies. Mediation analyses in clustered settings are planned
#' against a number of clusters, and the cluster counts that are actually
#' attainable in a real study are often much smaller than a study's nominal
#' sample. This reference table records both, so that planning can start from
#' documented figures rather than round numbers.
#'
#' @format A data frame with 3 rows and 10 columns:
#' \describe{
#'   \item{study}{Character. Study or data-source name.}
#'   \item{cohort}{Character. Cohort label, or `NA` when the source is not
#'     organized by cohort.}
#'   \item{clusters_sampled}{Integer. Number of clusters (schools) in the drawn
#'     sample, or `NA` when the source does not document one. For
#'     ECLS-K:2011 this is the original school sample drawn *before*
#'     substitution; see `sampling_note`.}
#'   \item{clusters_participating}{Integer. Number of clusters that actually
#'     participated.}
#'   \item{units_total}{Integer. Total number of level-1 units (children or
#'     students).}
#'   \item{mean_units_per_cluster}{Numeric. `units_total /
#'     clusters_participating`. Derived here, not taken from the source; note
#'     that the denominator is the *participating* cluster count, not the
#'     sampled count.}
#'   \item{cluster_response_rate}{Numeric. Cluster-level response rate as
#'     reported by the source, or `NA` when none is given. The ECLS-K:1998-99
#'     value is the *weighted* school response rate.}
#'   \item{restricted_use}{Logical. `TRUE` when a restricted-use data licence
#'     is required to obtain cluster (school) identifiers, which a multilevel
#'     analysis needs.}
#'   \item{sampling_note}{Character. Source-specific qualifications: which
#'     counts are approximate, which alternative counts the source reports,
#'     and why a field is `NA`.}
#'   \item{source}{Character. Citation, with a URL where the source is
#'     available online.}
#' }
#'
#' @details
#' NCES names the second cohort's study "ECLS-K:2011"; its `cohort` value here
#' is `"2010-11"`, the kindergarten class year. The two labels refer to the
#' same cohort.
#'
#' Two points about how the numbers should be read.
#'
#' *Participating versus sampled.* For both ECLS-K cohorts the participating
#' cluster count is well below the sampled count, and it is the participating
#' count that describes the data a user would actually analyze. The NCES
#' Handbook of Survey Methods summary page for ECLS-K:2011 implies roughly
#' 1,310 schools; the base-year user's manual (NCES 2015-074) gives roughly
#' 970 participating schools. The two figures come from different documents
#' and are not interchangeable. This table carries the user's-manual figure.
#'
#' *Approximate counts.* The ECLS-K:2011 participating-school and child counts
#' are reported as approximate in the user's manual and are stored here as
#' plain numbers; treat them as approximate.
#'
#' @source
#' ECLS-K:1998-99: NCES Handbook of Survey Methods,
#' \url{https://nces.ed.gov/statprog/handbook/ecls_k_dataquality.asp}
#'
#' ECLS-K:2011: User's Manual for the ECLS-K:2011 Kindergarten Data File
#' (NCES 2015-074), \url{https://nces.ed.gov/pubs2015/2015074.pdf}
#'
#' `mediation::student`: the \pkg{mediation} R package.
#'
#' @examples
#' multilevel_designs[, c("study", "cohort", "clusters_participating")]
#'
#' # Sources obtainable without a restricted-use licence
#' subset(multilevel_designs, !restricted_use)
"multilevel_designs"
