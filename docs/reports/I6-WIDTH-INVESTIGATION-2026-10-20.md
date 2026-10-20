# I6 finite-K width investigation (2026-10-20)

**Disposition:** The retained I6 data support finite-`K` Monte Carlo noise as
the explanation for candidate `d`'s s6/`N = 1000` interval-width excess.
This is a diagnostic finding, **not** a replacement for the frozen I6 gate.
I6 remains **FAIL**; I7, package-method claims, release, and pmed submission
remain on hold. The earlier [I6 gate report](I6-GATE-2026-10-19.md) and
[s6 diagnostic](I6-S6-N1000-DIAGNOSTIC-2026-10-19.md) remain authoritative.

## Retained input and method

The [plan](../../tasks/PLAN-i6-width-investigation-2026-10-20.md) uses no new
simulated data. The 300 retained Hopper job `4341753` CSV parts were copied
read-only from `/users/dtofighi/probmed_i6_20261019_full/parts/`. The local
and remote SHA-256 hash-list digests both equal
`658c2d583bd795cfdf9743cbaf8dbb37a1b8336845a9743a4a4cab17d4f59297f`.
The committed I6 runner hash is
`ac10f1a7440235aca3e834e8825c1507312c5d4219045549d895d81fe9e99adf`.
The diagnostic checked 60,000 unique keys, all 2,000 replication IDs in each
of 30 cells, 2,000 attempted bootstrap resamples per row, finite ordered
Gaussian intervals, and zero failed Gaussian `d` resamples.
All 20 recalculated candidate/exact coverages and width ratios match the
frozen I6 decision CSV exactly at stored precision.

Run `Rscript inst/sim/i6_width_investigation.R <parts_dir> <output.csv>` from
the repository root. The script pairs `d` and exact `e` within each outer
replication, resamples those paired rows 5,000 times with seed `20261020` for
a descriptive 95% interval on the ratio of mean widths, and writes the
[20-cell results](I6-WIDTH-INVESTIGATION-2026-10-20.csv). It predicts finite-`K`
broadening by adding `mean[p(1-p)] / 200000` to an exact-bootstrap variance
proxy inferred from the mean exact percentile width; `p` is each saved exact
bootstrap median. This normal-spread approximation is not a direct proof of
the endpoint distribution.

## Findings

At s6/`N = 1000`, the observed mean-width ratio is **1.029588** (descriptive
paired 95% CI **[1.028847, 1.030327]**), versus a predicted **1.029688**.
The lower endpoint moves down by `0.00026571` on average and the upper moves
up by `0.00026252`; mean midpoint shift is `-0.00000159`. Candidate `d` is
wider in 1,936 of 2,000 paired rows. There are 18 `d`-only coverage rows and
zero `e`-only rows, reproducing coverage 0.9515 versus 0.9425. The frozen
coverage-equivalence CI and width gate both still fail.

All 20 Gaussian cells have observed width ratios above 1. The finite-`K`
prediction tracks their increase as exact intervals narrow: the largest
absolute observed-minus-predicted ratio is `0.000519`. The prediction falls
just outside the descriptive ratio CI in 3 of 20 cells, so it should not be
treated as an exact interval model or as a calibrated inferential test.
Neither the near-symmetric endpoint movement nor the complete paired data
suggests a collation failure or a one-sided endpoint shift. The saved rows do
not contain every bootstrap draw, so they cannot exclude every within-row
mechanism.

**Next scientific decision:** Keep I6 FAIL. A Gaussian exact-calculation
method could be evaluated as a *new* candidate under a reviewed implementation
and prespecified validation plan; alternatively, a stochastic `K` proposal
would need its own resources and full-cell gate. These post-result diagnostics
must not be used to tune `K` and relabel I6 as a pass. Binary s3/s4 remain
report-only and unaddressed here.
