# Literature review: the retained exact-Gaussian coverage gate (2026-10-20)

**Verdict.** The original rule is a statistically coherent, conservative *equivalence-style* test of true coverage in six specified Gaussian null cells, provided the 6,000 coverage indicators in each cell are independent Bernoulli outcomes from the fixed data-generating procedure and all attempted runs are counted. Statistical-methods authors recommend prespecified protocols, Monte Carlo uncertainty, equivalence bounds justified by their application, and joint examination of coverage with interval width and failures. **The reviewed sources do not recommend 93%–97% as a universal tolerance for nominal 95% coverage.** Its substantive acceptability remains an author judgment. Retain the frozen gate for the completed run, but describe its acceptance only as “the six specified cells met the prespecified 93%–97% coverage criterion,” not as general calibration or proof that 93% coverage is harmless.

## Exact question and scope

The [retained plan](../PROPOSAL-exact-gaussian-followup-2026-10-19.md) requires, in each s6/s7 null cell at N = 250, 500, and 1,000, a *two-sided 90% Clopper–Pearson interval* for the coverage proportion to lie strictly inside (0.93, 0.97). Each cell has 6,000 outer replications and 2,000 nonparametric bootstrap draws per outer data set. This is a test about the coverage probability of a **95% percentile-bootstrap interval**, not a 90% interval for the mediation estimand. The separate probmed I6 candidate-method comparison has different frozen paired coverage/width gates and remains **FAIL**; this review cannot reverse it or authorize I7, package changes, or pmed submission.

## What the statistical literature supports

1. **The inferential form is valid, conditional on the margin.** Two one-sided tests (TOST) at 5% correspond to placing the relevant two-sided 90% interval inside prespecified equivalence bounds; Lakens, Scheel, and Isager explain this rule and emphasize substantive choice of bounds. Schuirmann developed its classical form for bioequivalence. Applying that logic to binomial coverage uses exact binomial bounds rather than their normal-mean intervals. The analogy supports the *test structure*, not the 93%–97% numbers. [Lakens et al. (2018)](https://journals.sagepub.com/doi/10.1177/2515245918770963); [Schuirmann (1987)](https://community.jmp.com/kvoqx44227/attachments/kvoqx44227/discussions/96942/1/TOSTarticle.pdf).
2. **All-six pass is a conjunction, not six discoveries.** For the global claim that *every* cell is inside its band, the null says *at least one* cell is outside. If each component test has size at most 0.05, the probability of falsely accepting the global claim is at most 0.05, regardless of dependence among cells: an all-six pass necessarily includes a false pass for any truly out-of-band cell. Intersection–union testing formalizes this. It does **not** mean the six reported 90% cell intervals are a simultaneous 95% confidence set, nor license individual post hoc “all cells certified” statements outside the prespecified global decision. [Berger and Hsu (1996)](https://www.public.asu.edu/~rlberge1/papers/statsci96.pdf).
3. **Simulation design and interpretation need more than a binary pass.** Burton et al. recommend a detailed simulation protocol; Morris, White, and Crowther recommend defining aims, mechanisms, estimands, methods, and performance measures in advance, choosing repetitions via Monte Carlo SE, and reporting uncertainty. Morris et al. distinguish undercoverage from overcoverage and note that coverage can conceal bias/width problems. Thus width, point-estimate behavior, failures, and out-of-gate cells need separate reporting and scrutiny, even when they do not change this frozen acceptance decision. [Burton et al. (2006)](https://onlinelibrary.wiley.com/doi/10.1002/sim.2673); [Morris et al. (2019)](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8086).
4. **The margin is a value judgment.** The reviewed methodological recommendations require bounds to reflect a meaningful tolerance; they do not supply a universal 93%–97% band for simulated 95% intervals. At the lower boundary, 7% of intervals miss versus 5% nominally—a 40% relative increase in misses. A symmetric ±2-point coverage band is therefore not a claim that under- and overcoverage have equal consequences. Its prior specification protects against tuning to this result, but does not itself prove substantive adequacy. [Lakens et al. (2018)](https://journals.sagepub.com/doi/10.1177/2515245918770963); [Morris et al. (2019)](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8086).

## Reproducible audit of the rule, not of the Hopper data

Let K be the number of covered outer intervals, n = 6,000, and p the fixed cell's true coverage. For 0 < K < n, the equal-tailed 90% exact limits are `qbeta(.05, K, n-K+1)` and `qbeta(.95, K+1, n-K)`. Strict containment in (0.93, 0.97) gives **5,613 ≤ K ≤ 5,797**, or observed coverage **93.55%–96.62%**. At the lower passing count, the lower limit is 0.9300417; at the upper passing count, the upper limit is 0.9699203. This was calculated with base R, independently of the run summary:

```r
n <- 6000L
k <- 0:n
lo <- qbeta(.05, k, n - k + 1)
hi <- qbeta(.95, k + 1, n - k)
range(k[lo > .93 & hi < .97])  # 5613 5797
```

Assuming independent outer draws *within a cell*, the exact one-cell pass probability is `pbinom(5797, 6000, p) - pbinom(5612, 6000, p)`: 0.0487 at p = 0.93, 0.9313 at p = 0.94, effectively 1 at p = 0.95, and 0.0463 at p = 0.97. The 95%-coverage Monte Carlo SE is `sqrt(.95*.05/6000) = 0.002814`, or **0.281 percentage points**. Under an *additional* independence assumption across six cells, all-six pass probability at p = 0.94 in every cell is about **0.652**; it may differ with shared random streams. Thus a nonpass near a boundary need not mean the procedure is bad, and a pass does not imply exactly 95% coverage. These probabilities audit the rule's operating characteristics, not the actual six coverage counts.

The pmed source-verification report at `/Users/davood/projects/research/pmed/docs/reports/EDIT-REPORT-pmed-exact-gaussian-2026-10-19.md` reports a frozen six-cell PASS, observed coverage 94.38%–95.32%, and zero incomplete intervals/nonfinite bootstrap estimates in those cells. That report says its saved summary was independently reaggregated. **This literature review did not itself access or reconcile the Hopper manifest, raw chunks, six integer counts, or hashes.** The subsequent [independent source audit](AUDIT-exact-gaussian-hopper-2026-10-20.md) completed that separate check.

## Method-by-setting evidence map

| Methodological question | One fixed Gaussian null cell | Six-cell global gate | Other models / cells |
| --- | --- | --- | --- |
| Prespecified protocol and Monte Carlo uncertainty | Burton; Morris et al.; applicable | Same, with all-six rule fixed in advance | Applicable as design principle, not validation |
| Equivalence-style CI containment | Schuirmann; Lakens et al.; exact-binomial adaptation derived here | Berger–Hsu intersection–union argument | Requires a new model-specific claim and justified margin |
| Coverage versus width, bias, and failures | Morris et al.; report separately | Do not average away weak cells | Report independently; no inference from six-cell pass |
| Particular 93%–97% tolerance | **No universal recommendation identified** in reviewed sources | Author's frozen criterion, not consensus | Not transferable |

## Claim-to-source traceability and limits

| Claim used here | Primary source or calculation | Limit |
| --- | --- | --- |
| Detailed, prospective simulation protocol | [Burton et al. (2006)](https://onlinelibrary.wiley.com/doi/10.1002/sim.2673); [Morris et al. (2019)](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8086) | Neither endorses this numerical band. |
| 90% interval inside bounds at one-sided α = .05 | [Schuirmann (1987)](https://community.jmp.com/kvoqx44227/attachments/kvoqx44227/discussions/96942/1/TOSTarticle.pdf); [Lakens et al. (2018)](https://journals.sagepub.com/doi/10.1177/2515245918770963) | Classical examples involve means/effects; binomial adaptation is our derivation. |
| All-component pass controls one global equivalence claim | [Berger and Hsu (1996)](https://www.public.asu.edu/~rlberge1/papers/statsci96.pdf) | Not simultaneous marginal CI coverage or broader model validity. |
| Count thresholds, Monte Carlo SE, pass probabilities | Base-R equations above | Assumes a fixed binomial coverage probability and within-cell independent attempts. |
| Coverage can mask width or bias | [Morris et al. (2019)](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8086) | Descriptive width review is not a retroactive gate. |

**Recommended decision.** Keep the original gate for this already completed experiment; do not substitute the later 94%-minimum proposal post hoc. After source reconciliation, report per-cell counts, coverage, exact limits, MCSE, widths, failures, and out-of-gate cells. If any saved-data integrity check or cell limit fails, withhold the scoped acceptance claim and investigate. A future independent study could choose a differently justified tolerance *prospectively*, but would be a new question, not a reinterpretation of this run.
