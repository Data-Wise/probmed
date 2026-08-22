# Literature review: the inference citations behind the gauge CI argument

**Created:** 2026-08-21
**Scope:** the four references load-bearing in `ward_residual()`'s CI machinery and in
the gauge-pmed manuscript's inference sections
**Companion:** `SPEC-2026-08-21-bootstrap-calibration.md` (the empirical side)

---

## Why this exists

Every prior characterization of these four papers in this project — in the manuscript
supplement, in issues #11 and #32, and in `R/gauge-pmed.R`'s own comments — traced
back to one-line summaries, not to the papers. None of the four is in the Zotero
library. This review reads them.

**Method.** Metadata verified against Crossref. Full texts retrieved and read.
Provenance is labeled throughout, because the point of the exercise is to stop
treating abstract-level knowledge as if it were full-text knowledge.

**What was NOT read, and matters:**

- **DiCiccio & Efron:** the published *discussion* (pp. 213-228) was not retrieved —
  only the main article, pp. 189-212. A percentile critique by the discussants would
  live exactly there.
- **Zhan:** read in the author's accepted manuscript (faculty page, dated 2026-05-18),
  not the Taylor & Francis version of record. The Online Supplementary Appendix A-F,
  containing all proofs and the many-instrument/wild-bootstrap simulations, was not
  retrieved.
- **Lin & Han:** read at arXiv:2604.17239v1; the SSRN copy returns 403, so the two
  versions could not be confirmed identical.

### Verification pass against the Zotero PDFs (2026-08-21, later same day)

PDFs were subsequently added to Zotero, allowing the agent-reported findings to be
checked against local copies. Results:

| Claim | Status |
|---|---|
| Lin & Han: "Without refitting the nuisance estimator" | **CONFIRMED** in the SSRN copy (Zotero `TKPBYGYN`, 30pp). "refit" occurs **exactly once** in the whole document; no simulation section. arXiv and SSRN agree on the decisive point — **gap 4 closed.** |
| DiCiccio & Efron: no "percentile method" | **CONFIRMED** — "percentile method/interval/bootstrap" = **0 occurrences** in the JSTOR PDF (`MEIT8R7K`). "necessary and sufficient" present as quoted. |
| Zhan: KS-distance mechanism, not interval comparison | **CONFIRMED** (`XPEQV45A`, 18pp): "instruments are deemed weak iff the Kolmogorov-Smirnov (KS) distance below exceeds..." |
| Zhan: anti-pretest footnote removed in publication | **CONFIRMED** — `pre-test`, `pretest`, `Guggenberger` all return **0 hits** in the published version. |

Two gaps narrowed but still open:

- **The D&E discussion is partially present.** The PDF runs to p. 213 and includes the
  opening page of **"Comment" by Peter Hall and Michael A. Martin** — but truncates
  after ~24 lines, before any substance. Notably the discussant is Hall, author of the
  Hall (1988) Edgeworth results that *both* DiCiccio & Efron and Owen (2025) build on.
  No "percentile" hits in the fragment available. Worth obtaining in full.
- **Zhan's Supplementary Appendix A-F is referenced but not included** ("Proof. See
  Supplementary Appendix A/B/C/D"; simulation evidence in E and F). All proofs remain
  unread. Available at the article DOI.

---

## 1. Zhan (2026), *Detecting weak identification by bootstrap*

Econometric Reviews **45**(7), 912-928. DOI 10.1080/07474938.2026.2670373.

### What the manuscript claims

Supplement, `gauge-pmed.qmd:422-424`: weak identification is *"signalled by divergence
between the Wald and percentile intervals [@zhan2026]."* Listed as validation gate
(ii). Implemented as probmed's A2 gate (`weak_id`): percentile-CI width over Wald-CI
width, flagged at >= 3x.

### What the paper actually does

**The diagnostic is not a comparison of two intervals — not by width, not by
endpoints, not by overlap.** No interval-to-interval comparison appears in the paper.

It compares **the bootstrap distribution of the standardized TSLS estimator against
the standard normal**, by Kolmogorov-Smirnov distance. Definition 1 *defines* weak
instruments as `KS = sup_c |P(sqrt(n)(theta-hat - theta)/sigma <= c) - Phi(c)| > 5%`.
The feasible version (eq. 11) uses a residual bootstrap with B = 10,000, centered at
`theta-hat` and scaled by the **conventional** `sigma-hat`; inference on KS comes from
a double bootstrap (§3.3); the test is `H0: KS <= 5%` (eq. 13).

The paper is explicit that bootstrap *failure* is the signal: "when the bootstrap
fails, the bootstrap distribution of the standardized TSLS estimator differs
substantially from the standard normal distribution, which signals the existence of
weak instruments."

### Three further scope problems

1. **Setting.** Classical linear IV/TSLS, i.i.d., homoscedastic, fixed `k`. Not GMM,
   not nonlinear, not semiparametric, and explicitly not cross-fit/DML. GMM appears
   once, in a literature footnote, as an *open demand* the paper does not meet.
2. **KS is confounded with endogeneity by construction.** Corollary 1's leading term
   is `rho*c^2*phi(c)/sqrt(mu^2)` — driven by the endogeneity correlation `rho` as
   well as the concentration parameter. §4.1.3: the F-test and the bootstrap test "do
   not simply substitute for each other nor are they directly comparable." Table 1:
   at `mu^2 = 5`, KS = 0.111 at `rho = 0.99` vs 0.079 at `rho = 0.50`. **KS is not a
   monotone function of identification strength, and Zhan says so.**
3. **The theory is a strong-identification theory.** Theorems 1-3 hold under
   Assumption 1 (Strong Instrument Asymptotics); the published version adds that
   "under weak instruments, the remainder terms of the Edgeworth expansion matter, so
   this paper only adopts the Edgeworth expansion under the null of strong
   instruments." Power against weak identification rests on a concentration-parameter
   *pattern* result plus Monte Carlo — **not** a consistency-of-detection theorem.

Ratios, Fieller, and cross-fitting appear nowhere.

### Verdict

**Misattributed on both mechanism and setting.** The manuscript attributes to Zhan a
Wald-vs-percentile width comparison in a DML ratio setting; Zhan proposes a
KS-distance-to-normal test in linear IV. Saying Zhan "proves" a diagnostic detects
weak identification would additionally overclaim.

### Consequence for our own validation result

This **reframes** `.STATUS` item 5's finding rather than contradicting it. The 48k
grid found A2 (~97% width-explained) carries little information beyond interval
narrowness. That is a fair verdict on **the gate as implemented** — but it was never
a test of Zhan's procedure, because the implemented gate is not Zhan's procedure. The
honest statement is now stronger and simpler: *the A2 gate is not Zhan's method, and
on its own terms it is a narrowness proxy.* Whether Zhan's actual KS test would work
here is untested and, given the setting mismatch, not obviously worth testing.

---

## 2. DiCiccio & Efron (1996), *Bootstrap confidence intervals*

Statistical Science **11**(3), 189-228. DOI 10.1214/ss/1032280214.

### What the manuscript claims

`gauge-pmed.qmd:415-419`: the percentile interval is "the appropriate construction for
a skewed ratio, for which a wider *symmetric* SE still mis-covers [@dicicco1996]."

### What the paper actually says

**The paper does not present a "percentile method."** The strings "percentile method,"
"percentile interval," and "percentile bootstrap" occur **zero times** in pp. 189-212.
The abstract names its four procedures: "BC_a, bootstrap-t, ABC and calibration."

The percentile interval appears exactly once — as the **degenerate case of BC_a** when
both corrections are switched off (p. 193, after eq. 2.3): "If a and z_0 are zero,
then theta-hat_BC[alpha] = G-hat^{-1}(alpha), the 100*alpha-th percentile of the
bootstrap replications."

The architecture is a ladder: **standard subset percentile subset BC_a**. And the
paper states the sufficiency claim as near-theorem (p. 194): the three extensions
(bias, nonconstant standard error, normalizing transformation) "are **necessary and
sufficient** to give second-order accuracy." The percentile interval keeps one of the
three — transformation invariance — and discards precisely the skewness correction.

### The clause that IS supported, and its better citation

The manuscript's *premise* is strongly supported, and stated in the manuscript's own
terms — that width cannot repair shape (p. 190):

> "**The standard intervals always have shape equal to 1.00. It is in this way that
> they err most seriously.** ... This kind of error is automatically identified and
> corrected by all the bootstrap confidence interval methods."

Cite **p. 190's shape decomposition** for that premise, not Table 2 alone.

### They have a worked ratio example — and its lesson

pp. 198, the cell data (n = 1,843): the parameter of interest is a **ratio of two
estimated probabilities**, `theta = pi_15/pi_51`. The result worth knowing:

> "Changing parameters from theta = pi_15/pi_51 to phi = log(theta) changes
> (a-hat, z_0-hat, c_q-hat) from (-0.006, -0.025, 0.105) to (-0.006, -0.025, 0.025)...
> **The standard intervals are nearly correct on the phi scale.** The ABC and BC_a
> methods automate this kind of data-analytic trick."

So on their own ratio example: the skewness breaking the symmetric interval is largely
a **scale artifact**, repairable by reparameterization. (Note this does not transfer
directly to `W = R/OE`, which takes negative values, so `log` is unavailable — but the
principle, reparameterize away from the singularity, does.)

Also worth noting: the paper's opening catalogue of the few *exact* intervals names
"the ratio of normal means" (p. 189) — the Fieller problem, named obliquely and not
pursued.

### What it does not say

No Fieller. No near-zero denominators. **No distinction between heavy tails and
skewness** — §8 assumes heavy tails away ("the fourth- and higher-order cumulants are
of order O(n^-1) or smaller," p. 203). The three "long-tailed" mentions are
descriptive and one-sided, and in each the point is that BC_a correctly extends
*further* into the tail — the opposite of a too-wide warning.

On BC_a stability the paper is *reassuring*, not cautionary: `a-hat`, `z_0-hat` err by
`O_p(n^-1)` and (p. 205, eq. 8.7) "this change does not affect the second-order
correctness." The instability warnings in the paper attach to **bootstrap-t** ("can be
numerically unstable, resulting in very long confidence intervals... a particular
danger in nonparametric situations," p. 199), not to BC_a — an easy misattribution.

### Verdict

**Premise supported; prescription inverted.** DiCiccio & Efron argue that being
tail-aware means *estimating* `z_0` and `a` — not reading raw quantiles. The
manuscript cites them for the rung of the ladder they treat as deficient.

But the repair is not simply "switch to BC_a." §8's machinery presupposes asymptotic
normality, Cornish-Fisher expansions, and twice-differentiable functionals. A
near-zero denominator can produce a bootstrap distribution that is bimodal or has no
finite variance. **In that regime this paper justifies neither percentile nor BC_a.
It is silent, and silence is not endorsement.**

---

## 3. Lin & Han (2026), *Bootstrap consistency for general double/debiased ML estimators*

SSRN preprint, DOI 10.2139/ssrn.6695948; read at arXiv:2604.17239v1.

Full treatment in `SPEC-2026-08-21-bootstrap-calibration.md` §Phase 0b. Summary:

**The paper proves consistency for a bootstrap that does not refit nuisances** (§2.2,
eq. 2.3): "**Without refitting the nuisance estimator** eta-hat_{0,k}...". Only the
fold-restricted empirical measure is reweighted; the fold partition is fixed across
resamples. The word "refit" occurs once in the paper — in the sentence excluding it.
There is no simulation section.

`ward_residual()` resamples rows and refits the entire cross-fit, all nuisances, every
resample. **Different estimator; Theorem 3.2 does not cover it.**

Ratios are never discussed. The manuscript's "requires `OE` bounded away from 0" is
recoverable by embedding the ratio as `psi = phi_R - theta*phi_OE`, giving
`J_0 = -OE`, so A3.2(iii)'s singular-value bound becomes `c_0 <= |OE| <= c_1` — a
derivation, not the paper's text, and one that also bounds `|OE|` **above**, which no
surface currently mentions.

**Correct citation for the refit case:** Tang & Westling, arXiv:2404.03064, which
permits a bootstrap nuisance refit (condition B2) and warns that refitting analogously
can fail when the learner is sensitive to the tied/replicated observations Efron's
bootstrap produces. Caveat: their framework uses Donsker-type conditions, which
cross-fitting exists to avoid. *(Reported from full text by the reading agent; verify
directly before citing.)*

### Verdict

**Miscitation of the central object** — the one finding here that most clearly needs
fixing before submission, since it is not a matter of emphasis but of which theorem
covers which procedure.

---

## 4. Owen (2025), *Better bootstrap t confidence intervals for the mean*

arXiv:2508.10083 [math.ST], v2 2025-08-20. Stanford. Preprint, no journal venue —
which is very likely why issue #11 could only cite it informally.

**Identified with high confidence.** A second candidate (arXiv:2501.07645, "Coverage
errors for Student's t confidence intervals...") is ruled out: it is cited *inside*
2508.10083 and mentions neither BC_a nor percentile.

### What issue #11 claims

"literature split on percentile-vs-BCa (DiCiccio & Efron 1996; Owen 2025 warns BCa
under-covers)."

### What the paper says

The BC_a claim is **literally accurate**: abstract, "The BCa bootstrap produces
shorter intervals but **tends to severely under-cover the mean**." Discussion: "we can
rule out BCa... They all have below nominal coverage in **every example considered
here**."

But the scope conditions are restrictive and the issue drops all of them:
- **Estimand: the univariate mean only.** Every result.
- **Sample size 2 <= n <= 20.** Both BC_a and bootstrap-t are second-order accurate as
  `n -> inf`; the undercoverage is explicitly a **small-n** phenomenon.
- Worse under skewness/kurtosis: Exp(1) at n=10 gives BC_a coverage 0.877.

### The inversion

**Owen rates the percentile method *below* BC_a** (l.581): "This is known as the
percentile method. **It is not competitive with the bootstrap t method** when we want
an ACI for the mean. **It is not as well regarded as the BCa method** described next."
He never includes percentile in his own simulations, and his actual recommendation is
a *third* method — a bootstrap-t with Beta(1/2, 3/2) weights.

He also relays (L'Ecuyer et al. 2023, 2400+ RQMC tasks): "The percentile bootstrap
failed 1698 times, the bootstrap t failed 81 times, and the plain Student's t interval
failed only 3 times."

### Ratios and heavy tails

No ratio estimands, none. Heavy tails appear as heavy tails **of the data** (t_4,
lognormal), not of a ratio's bootstrap distribution. The Hall (1988) coefficients he
tabulates make BC_a ~3.6x more kurtosis-sensitive than percentile — but that expansion
assumes finite **eighth** moments plus a Cramer condition, which a ratio with a
near-zero denominator may flatly violate.

The most transferable piece is structural (l.2188): "the ACIs from BCa are nested
inside the convex hull [x_(1), x_(n)] of the data and that already limits their
coverage for small n."

### Verdict

**Accurate about BC_a; misrepresents the paper's position in the debate it is cited
into.** Issue #32's audit used Owen to argue against BC_a while proposing percentile —
against the method Owen prefers of the two. The claim should be scoped to "for the
univariate mean at small n" or dropped.

---

## Synthesis: the pattern

Four references, four different problems, all pointing one way.

| Reference | Cited for | What it actually supports |
|---|---|---|
| Zhan 2026 | Wald-vs-percentile divergence detects weak ID | A KS-to-normal test in linear IV; KS confounded with endogeneity; strong-ID theory |
| DiCiccio & Efron 1996 | Use the percentile interval for a skewed ratio | That symmetric intervals mis-cover skewed estimands. Percentile is their *deficient* rung |
| Lin & Han 2026 | Bootstrap consistency for our refit-per-resample scheme | Consistency for a **no-refit** reweighting scheme |
| Owen 2025 (issue #32) | BC_a under-covers, so prefer percentile | BC_a under-covers **the mean at n <= 20**; he rates percentile lower still |

**The percentile bootstrap arm has no citation support in this manuscript.** Its three
supporting references either argue for something else (DiCiccio & Efron), cover a
different procedure (Lin & Han), or address a different question in a different model
class (Zhan).

This converges with the independent empirical finding (SPEC Phase 0): the percentile
arm covers at **1.00 in 8 of 8 simulation cells**. Theory and evidence agree on which
component is the weak link.

**And no citation in this set covers the regime that actually matters here.** Every
one of the four either assumes away, or never considers, a ratio whose denominator can
approach zero. DiCiccio & Efron assume `O(n^-1)` fourth cumulants; Owen's Hall
expansion needs eighth moments and a Cramer condition; Lin & Han's A3.2(iii) requires
the Jacobian bounded away from singularity — which *is* the condition being violated.
The near-null regime is unsupported by all four.

## What this implies

1. **Fieller gains support from an unexpected direction.** DiCiccio & Efron's own
   catalogue of exact intervals names "the ratio of normal means" (p. 189). Combined
   with the exact identity established separately — A1's gate `oe_snr >= 2` **is** the
   Fieller boundedness criterion, since `seOE^2 = VOE` exactly, so
   `fa > 0 <=> oe_snr > 1.959964` — the case for reporting Fieller for `W` rests on
   the ratio literature rather than on the bootstrap literature that does not cover
   this regime.
2. **The no-refit multiplier bootstrap is now doubly motivated** — it is both the
   cheap discriminating experiment (SPEC Phase 0b) and the only bootstrap in this
   citation set actually proven for a cross-fit DML functional.
3. **Reparameterization deserves consideration** (DiCiccio & Efron p. 198). `log` is
   unavailable for a sign-changing `W`, but the principle stands.
4. **BC_a (#32) is not rehabilitated by this review.** DiCiccio & Efron justify it
   over percentile *in their regime*; neither is justified in ours. It remains the
   wrong next move, now for a documented reason rather than an assumed one.

## Corrections owed to project surfaces

- **Manuscript** (cross-repo, needs go-ahead): Zhan mechanism + setting; Lin & Han
  refit mismatch; DiCiccio & Efron clause re-scoped to the premise.
- **`R/gauge-pmed.R`** in-code comments and roxygen: the A2 gate is described as
  implementing "Zhan 2026." It does not. The warning text `weak_id` emits also names
  Zhan.
- **Issue #11**: the A1/A2 hierarchy note is correct as written, but its framing
  ("Zhan's premise did not survive") should become "the implemented gate is not
  Zhan's method."
- **Issue #32**: the Owen characterization needs scoping.
- **`.bib`**: `lin2026` authors and type; `zhan2026` volume/issue/pages.

## Unverified, and worth closing

Status after the Zotero verification pass (see §Why this exists):

1. **DiCiccio & Efron's discussion, pp. 213-228** — PARTIALLY OPEN. The local PDF
   reaches only the first page of Hall & Martin's Comment. Hall is the author of the
   Edgeworth results both this paper and Owen (2025) build on, so his comment is the
   single most valuable unread item in this set. Obtain from JSTOR/Project Euclid.
2. **Zhan's Supplementary Appendix A-F** — OPEN. Referenced throughout ("Proof. See
   Supplementary Appendix A"), not bundled with the article PDF. Contains every proof
   plus the many-instrument and wild-bootstrap simulations. (The related question of
   whether the published version restores the 2017 draft's anti-pretest footnote is
   now ANSWERED: it does not — zero hits for `pre-test`/`pretest`/`Guggenberger`.)
3. **Tang & Westling** (arXiv:2404.03064) — OPEN, and now the main gap for the
   constructive path. Metadata verified (Zhou Tang & Ted Westling; v1 2024-04-03,
   v2 2024-04-18; no journal reference), but the content is known only via the reading
   agent. It is the proposed replacement citation for the refit-per-resample
   bootstrap, so it should be read directly before anything is written on its
   authority. Not currently in Zotero.
4. ~~Whether the arXiv and SSRN versions of Lin & Han are identical~~ — **CLOSED.**
   The decisive no-refit sentence, the single occurrence of "refit", and the absence
   of a simulation section all reproduce in the SSRN copy.
