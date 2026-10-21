<!-- Subagent literature report (savant literature-gap-finder method), 2026-10-21. Model output, not verified by the author. OPENED = the subagent read the text; SEARCH-ONLY = existence confirmed from snippets only. Treat as data. -->

# Literature check: cusp functional P = sqrt|d1| / (sqrt|d1| + sqrt|d2|)

Status key. OPENED = I read the text (pdftotext of the PDF) and quote it. SEARCH-ONLY = confirmed to exist and summarized from search-result text, not opened; no quotes attributed. DERIVED = my own elementary argument, not a citation; verify before use.

## (i) Claim table

| # | Claim | Supported? | Source | Quote or [CITE-NEEDED] |
|---|---|---|---|---|
| 1a | Two-point (Le Cam) lower bound: risk floor controlled by total variation between two nearby distributions | Standard method; Tsybakov book NOT opened | Tsybakov (2009), Intro to Nonparametric Estimation, Springer, Ch. 2 | [CITE-NEEDED: theorem number and TV form not opened]. Opened secondary: Cai & Low (2011, Ann. Statist. 39(2):1012-1041): "Le Cam's method is a well-known approach based on this idea." and "depend on the total variation distance between the convex hulls" (Le Cam 1973, 1986, as cited there). |
| 1b | Donoho-Liu: minimax bounds for functionals via testing | Supported (as cited by Cai-Low) | Cai & Low 2011, Sec. 7 | "bounds on minimax estimation can be based on testing two composite hypotheses" (attributed to Donoho and Liu 1991). Donoho-Liu itself NOT opened. |
| 1c | Bias/SD two-point inequality, (E1 T - E0 T)^2 <= chi2(P1,P0) Var0(T) | Standard (Hammersley/Chapman-Robbins) | not opened | [CITE-NEEDED]; derivable in 3 lines (Cauchy-Schwarz on the likelihood ratio). |
| 2 | Hodges superefficiency; uniform/local-minimax frame is correct; lone-cell minimax is trivially 0 | Not opened. Logic is elementary (constant estimator has zero risk at its own point). | Hodges (1951 unpublished, per Le Cam 1953); Le Cam (1953); Hajek (1972) | [CITE-NEEDED] for all three. Opened related: Hirano-Porter state the target "regular estimators" notion and that nonregularity blocks convolution-theorem optimality (see row 3). |
| 3a | No locally asymptotically unbiased and no regular estimators if the functional is not continuously differentiable (single point suffices) | Supported, OPENED (2009 MPRA draft; published as Econometrica 80(4):1769-1790, 2012) | Hirano & Porter 2012 | Abstract: "there exist no locally asymptotically unbiased estimators and no regular estimators." Also: "nondifferentiability of κ̇(h) at a single point h0" suffices; Theorem 2: "Tn is not locally asymptotically unbiased and is not regular." Their remark: "bias correction procedures cannot completely eliminate local bias, and reducing bias too much will eventually cause the variance of the procedure to diverge." Note the paper is about unbiasedness/regularity; it states no rate. Published-version mention of local asymptotic minimaxity is SEARCH-ONLY. |
| 3b | Bias can be made arbitrarily small only at the cost of variance blowing up | Supported via Hirano-Porter's citation (OPENED) | Doss & Sethuraman 1989, Ann. Statist. 17(1):440-442 | H-P: "such a sequence must have variance increasing to infinity at every point in the parameter space." Condition is pointwise bias -> 0. Doss-Sethuraman itself SEARCH-ONLY (abstract: no unbiased estimator + biases -> 0 implies variances -> infinity). |
| 3c | Follow-up on the bias-variance tradeoff | SEARCH-ONLY | Jeesen Chen (2004), "Notes on the bias-variance trade-off phenomenon," IMS Lecture Notes Monograph, A Festschrift for Herman Rubin. Author and venue are from the search listing; the user's "Chen (2004)" is probably this. | [CITE-NEEDED: not opened]. Search summary: uses Hellinger distance to bound the mean of an estimator by its variance. |
| 3d | Informative unbiased estimators do not exist in singular problems (rate slower than n^-1/2) | SEARCH-ONLY | Liu & Brown 1993, Ann. Statist. 21(1):1-13 | [CITE-NEEDED: not opened]. Search abstract: informative unbiased estimators "do not exist in singular problems." Relevant because sqrt\|d\| has an n^-1/4 rate. |
| 3e | Pötscher-Schneider / thresholding non-uniformity | NOT confirmed | Pötscher & Leeb (2009), J. Multivariate Anal. 100:2065-2082 exists (hard-thresholding, LASSO, SCAD; "uniform convergence rate" derived; SEARCH-ONLY). | The "Pötscher and Schneider (2010)" item the project names was not located. [CITE-NEEDED]. Do not cite a 2010 Pötscher-Schneider thresholding paper until verified. |
| 3f | Dümbgen (1993) on nondifferentiable functionals and the bootstrap | Supported via Fang-Santos (OPENED) | Dümbgen 1993 (full title and venue not opened) | Fang-Santos: Dümbgen "first examined the validity of the bootstrap for estimating the asymptotic distribution of φ(θ̂n) under a potential lack of differentiability." |
| 4a | Standard bootstrap inconsistent iff the map is not (fully) differentiable (Gaussian limit) | Supported, OPENED | Fang & Santos 2019, Rev. Econ. Stud. 86(1):377-412 | Theorem 3.1: φ "is (fully) Hadamard differentiable ... if and only if" the bootstrap consistency display holds. Their Example 2.1 is exactly \|E[X]\|. |
| 4b | Alternatives that stay valid | Supported, OPENED | same | "alternative resampling schemes such as the m out of n bootstrap ... and subsampling can nonetheless remain valid"; they propose an estimated-directional-derivative scheme with local size control. |
| 4c | Numerical bootstrap | Partly. Fang-Santos cite "Hong and Li (2015)": the directional derivative "can often be obtained through numerical differentiation of φ at θ̂n." Final paper: Hong & Li, Ann. Statist. 48(1):397-412 (2020) (SEARCH-ONLY, venue/pages from search). | | Quote above is from the Fang-Santos text, attributing to the 2015 working version. Final paper's own text NOT opened. |
| 4d | Andrews (2000, Econometrica) boundary | NOT opened | Andrews 2000, Econometrica 68(2) ("Inconsistency of the bootstrap when a parameter is on the boundary of the parameter space") | [CITE-NEEDED]. Relevant analogy only; it is not about a sqrt\|.\| cusp. |
| 4e | Valid SEs near the cusp: delta-method SE is not valid at d=0 | Supported by Fang-Santos (delta method for directionally differentiable maps gives a non-Gaussian limit) | Fang & Santos | Opened Theorem 3.1 shows the naive bootstrap approximates the wrong law. They give no "SE" per se; the limit law is non-normal, so a single SE is not a sufficient summary. |
| 5a | Fieller sets can be intervals, complement of an interval, or the whole line | Supported, OPENED | von Luxburg & Franz 2009, Statistica Sinica 19(3):1095-1117 | "two disjoint unbounded intervals, or even of the whole real line." |
| 5b | Any method that cannot output unbounded sets has arbitrarily bad coverage (Gleser-Hwang) | Supported via von Luxburg-Franz (OPENED); Gleser-Hwang itself SEARCH-ONLY | Gleser & Hwang 1987, Ann. Statist. 15(4):1351-1362 | vL-F: "any method which is not able to generate such unbounded confidence limits for a ratio leads to arbitrary large deviations from the intended confidence level." Also: "bootstrapping ρ̂ directly is not an option." |
| 6a | Proportion mediated is unstable when the total effect is near zero | SEARCH-ONLY, informal | MacKinnon et al. (as cited in a 2026 preprint, arXiv 2609.06445); Siani et al. article; Franz (2007) cited in arXiv 1612.08699 | [CITE-NEEDED]: no primary mediation source opened. |
| 6b | Wasserstein-based mediation effects with a cusp at zero effects | NOT found. Related but not the estimand: arXiv 2208.13939 (density mediators, barycenters), arXiv 2603.15182 (sequential transport), arXiv 2602.19290 (Wasserstein distance as a distributional RD effect). | | [CITE-NEEDED]. No paper found stating non-differentiability at zero for a Wasserstein proportion. Treat the project's cusp result as novel unless a closer paper is found. |
| 7a | Minimax rate for a nonsmooth functional of a Gaussian mean | Supported but NOT the same problem | Cai & Low 2011 | OPENED: "the minimax risk equals β*² M² (log log n / log n)² asymptotically" for (1/n) Σ\|θ_i\|, Y ~ N(θ, I_n). This is the average over n coordinates (dimension = n, noise level 1), NOT a single \|θ\| with n replicates. The log rate does NOT transfer to your 1-D d_j problem. |
| 7b | Minimax rate n^-1/4 for sqrt\|d\| in a 1-D Gaussian shift | DERIVED, no source found | see Derivation below | [CITE-NEEDED]. No opened source states it. |
| 7c | Lepski-Nemirovski-Spokoiny (1999) L1-norm of a regression function | Existence confirmed via Cai-Low | | Cai-Low: LNS "attention was focused on estimating the L1 norm of a regression function." LNS itself NOT opened. |

## Derivation for the project (7b, 1a, 1c). Mine, not cited. Check before using.

Setup: d̂ ~ N(d, s²) with s = σ/√n; g(d) = sqrt|d|.
- Upper bound: |sqrt|a| - sqrt|b|| <= sqrt(|a-b|), so E|g(d̂) - g(d)| <= E sqrt|d̂-d| <= s^(1/2) for all d, and RMSE-type versions are of the same order. The plug-in sqrt|d̂| attains order n^-1/4 uniformly.
- Lower bound (Le Cam two-point): take d0 = 0 and d1 = c·s. Then |g(d1) - g(d0)| = c^(1/2) s^(1/2), and TV(N(0,s²), N(cs,s²)) = 2Φ(c/2) - 1 < 1, fixed. Hence inf over estimators of max over {d0,d1} of risk >= const · s^(1/2) = const · n^-1/4.
- Conclusion supported by the argument: sup risk over a 1/√n neighborhood of 0 is of order n^-1/4 for EVERY estimator. So "no threshold improves the n^-1/4 rate" is true for worst-case (local minimax) risk, and false pointwise at d = 0 (a constant or threshold estimator has risk 0 at exactly 0; this is Hodges). The project must say "local minimax," never "at the point 0."
- Bias/SD gate in two nearby cells: for any estimator T, (E1 T - E0 T)² <= χ²(P1,P0)·Var0(T). With |bias_j| <= 0.10·SD_j in both cells this forces |g(d1) - g(d0)| <= (0.2 + sqrt(χ²))·SD, so two distinguishable-from-nothing cells whose targets differ by much more than the SD cannot both pass. Caveats: (1) the gate is scale-free, so an estimator can pass by inflating SD; this is the Doss-Sethuraman mechanism and argues for gating RMSE, not bias/SD alone. (2) The statement is conditional on Δg/SD, not an unconditional impossibility. For P (not just g), the two-cell Δ must be computed on P itself.
- P is a smooth function of (g1, g2) away from g1 = g2 = 0 only (ratio of nonnegative terms; at d1 = d2 = 0 the denominator vanishes too). So the P-level cusp has a second singularity at the joint origin; the rate n^-1/4 for P near a single zero contrast d1 = 0 with d2 fixed away from 0 follows from the g-level rate by the delta method in (g1, g2). At d1 = d2 = 0 P is undefined and the gate must exclude it or define it.

## (ii) What the project may safely cite / must not claim

May safely cite (all opened):
- Hirano & Porter (2012): no locally asymptotically unbiased or regular estimators when the functional is not continuously differentiable; one bad point suffices; bias reduction "will eventually cause the variance ... to diverge" (this last phrasing is theirs, based on Doss-Sethuraman).
- Fang & Santos (2019): standard bootstrap consistent iff full Hadamard differentiability (Gaussian limit); m-out-of-n and subsampling can remain valid; directional-derivative-based alternative scheme.
- von Luxburg & Franz (2009): Fieller sets can be unbounded; Gleser-Hwang as cited there.
- Cai & Low (2011): for context, as the example that nonsmooth functionals (|θ|) can have slower-than-parametric minimax rates, but only for the n-dimensional ℓ1 average.

May cite with the "per secondary source" tag: Doss & Sethuraman (1989) (through Hirano-Porter); Gleser & Hwang (1987) (through von Luxburg-Franz); Dümbgen (1993) (through Fang-Santos); Hong & Li (through Fang-Santos).

Must not claim:
- That any opened source proves "no estimator can meet a bias/SD < .10 gate." No source does; it is a project-specific two-point argument and must be proved in the paper.
- That the minimax rate of sqrt|d| is n^-1/4 "by Cai and Low" or "by Hirano and Porter." Neither states it. Cai-Low's (log log n / log n)² is a different problem.
- That "no threshold improves n^-1/4" without "local minimax / sup over a neighborhood of 0."
- That Hirano-Porter gives a rate or a minimax lower bound; the opened text proves nonexistence of unbiased/regular estimators only.
- Any "Pötscher and Schneider (2010)" thresholding reference, until located. Do not cite Andrews (2000), Tsybakov theorem numbers, Hodges, Le Cam 1953, Hajek 1972, Chen 2004, Liu-Brown beyond "exists," or any Wasserstein-mediation cusp paper from this review.
- That Fieller intervals have exact coverage for P. Fieller-type exactness is for ratios of normal means; P is not a ratio of linear combinations of means.

## (iii) Gaps already settled by the literature, and open ones

Settled:
1. Non-existence of unbiased/regular estimators at a cusp (Hirano-Porter 2012), so a "zero bias" gate at exact zero is unachievable in the unbiased-estimator sense and bias reduction costs variance.
2. Standard bootstrap fails at the cusp; m-out-of-n, subsampling, and directional-derivative schemes are the established fixes (Fang-Santos 2019). A percentile bootstrap SE at d near 0 is a known-bad choice; the project's Fieller-type inversion is consistent with that literature.
3. Unbounded Fieller sets are necessary for coverage (Gleser-Hwang as cited by von Luxburg-Franz). If the project's union-of-lines interval is always bounded, flag it.

Open (no opened source settles it):
- The n^-1/4 local-minimax rate for sqrt|d| (derived above in about 5 lines; needs a clean reference or a proof in the paper).
- Whether hard thresholding yields a strictly smaller constant, not rate, in the local-minimax sense.
- The weak-signal ladder risk gate calibrated to the lower bound (design step (b)); the constant in the two-point bound is explicit: best two-point constant comes from optimizing c in the TV bound.
- Any primary source on the proportion-mediated instability with a quotable statement, and any Wasserstein-mediation paper with this cusp.

## Not opened (could not access text)
Tsybakov (2009); Le Cam (1973, 1986, 1953); Donoho-Liu (1991); Hodges; Hajek; Doss-Sethuraman (Project Euclid returned an empty page, curl blocked by Incapsula); Chen (2004); Liu-Brown (empty page); Gleser-Hwang (empty page); Andrews (2000); Dümbgen (1993); Hong-Li (2020; Project Euclid empty); Pötscher-Leeb/Schneider; Lepski et al. (1999); MacKinnon et al. on proportion mediated; all Wasserstein-mediation items.
