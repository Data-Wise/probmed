# Design metadata for candidate multilevel studies

Publicly documented design metadata for a small set of clustered
(multilevel) studies. Mediation analyses in clustered settings are
planned against a number of clusters, and the cluster counts that are
actually attainable in a real study are often much smaller than a
study's nominal sample. This reference table records both, so that
planning can start from documented figures rather than round numbers.

## Usage

``` r
multilevel_designs
```

## Format

A data frame with 3 rows and 13 columns:

- study:

  Character. Study or data-source name.

- cohort:

  Character. Cohort label, or `NA` when the source is not organized by
  cohort.

- design_type:

  Character. `"probability sample"` for a designed sample with a
  documented frame, `"teaching data"` for a convenience data set shipped
  for instruction.

- clusters_sampled:

  Integer. Number of clusters (schools) in the original drawn sample,
  *before* any substitution, or `NA` when the source does not document
  one.

- clusters_eligible:

  Integer. Number of sampled clusters found eligible for data
  collection, or `NA` when the source does not document one.

- clusters_recruited:

  Integer. Total number of clusters approached, including any
  substitutes added after the original sample, or `NA` when the source
  does not document one.

- clusters_participating:

  Integer. Number of clusters that actually participated. This count may
  include substitute clusters, so it does not share a base with
  `clusters_sampled`.

- units_total:

  Integer. Total number of level-1 units (children or students).

- mean_units_per_cluster:

  Numeric. `units_total / clusters_participating`. Derived here, not
  taken from the source; note that the denominator is the
  *participating* cluster count, not the sampled count.

- cluster_response_rate:

  Numeric. Cluster-level response rate as reported by the source, or
  `NA` when none is given. The ECLS-K:1998-99 value is the *weighted*
  school response rate.

- restricted_use:

  Logical. `TRUE` when a restricted-use data licence is required to
  obtain cluster (school) identifiers, which a multilevel analysis
  needs.

- sampling_note:

  Character. Source-specific qualifications: which counts are
  approximate, which alternative counts the source reports, and why a
  field is `NA`.

- source:

  Character. Citation, with a URL where the source is available online.

## Source

ECLS-K:1998-99: NCES Handbook of Survey Methods,
<https://nces.ed.gov/statprog/handbook/ecls_k_dataquality.asp>

ECLS-K:2011: User's Manual for the ECLS-K:2011 Kindergarten Data File
(NCES 2015-074), <https://nces.ed.gov/pubs2015/2015074.pdf>

[`mediation::student`](https://rdrr.io/pkg/mediation/man/student.html):
the mediation R package.

## Details

NCES names the second cohort's study "ECLS-K:2011"; its `cohort` value
here is `"2010-11"`, the kindergarten class year. The two labels refer
to the same cohort.

Three points about how the numbers should be read.

*Recruitment chain.* `clusters_sampled`, `clusters_eligible`,
`clusters_recruited`, and `clusters_participating` are successive stages
with different bases. For ECLS-K:2011, 1,352 schools were sampled, 1,264
were eligible, 93 substitutes brought recruitment attempts to 1,446, and
roughly 970 participated. Because the participating count includes
substitutes, `clusters_participating / clusters_sampled` is *not* a
response rate; `cluster_response_rate` is `NA` on that row for this
reason, and only the rate a source reports directly is carried.

*Participating versus sampled.* For both ECLS-K cohorts the
participating cluster count is well below the sampled count, and it is
the participating count that describes the data a user would actually
analyze. The NCES Handbook of Survey Methods summary page for
ECLS-K:2011 implies roughly 1,310 schools; the base-year user's manual
(NCES 2015-074) gives roughly 970 participating schools. The two figures
come from different documents and are not interchangeable. This table
carries the user's-manual figure.

*Approximate counts.* The ECLS-K:2011 participating-school and child
counts are reported as approximate in the user's manual and are stored
here as plain numbers; treat them as approximate.

## Examples

``` r
multilevel_designs[, c("study", "cohort", "clusters_participating")]
#>                study  cohort clusters_participating
#> 1             ECLS-K 1998-99                    940
#> 2             ECLS-K 2010-11                    970
#> 3 mediation::student    <NA>                    568

# Sources obtainable without a restricted-use licence
subset(multilevel_designs, !restricted_use)
#>                study cohort   design_type clusters_sampled clusters_eligible
#> 3 mediation::student   <NA> teaching data               NA                NA
#>   clusters_recruited clusters_participating units_total mean_units_per_cluster
#> 3                 NA                    568        9679               17.04049
#>   cluster_response_rate restricted_use
#> 3                    NA          FALSE
#>                                                                                                                            sampling_note
#> 3 Public teaching data set shipped with the mediation R package. No sampled-school count or response rate is documented, so both are NA.
#>                                                         source
#> 3 mediation R package (`data(student, package = "mediation")`)
```
