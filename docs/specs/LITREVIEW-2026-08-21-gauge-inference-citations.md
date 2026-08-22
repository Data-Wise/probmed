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

**Retrieval history.** The first pass worked from whatever full texts could be reached
on the open web, and three gaps were flagged as unread: DiCiccio & Efron's published
discussion, Zhan's Online Supplementary Appendix, and confirmation that Lin & Han's
SSRN and arXiv versions agreed. **All three were subsequently closed** as the sources
were added to Zotero over the course of the session. Two of the closures changed
conclusions in this review rather than merely confirming them — see the Hall & Martin
correction and the rejoinder finding in §2. The sequence is left visible below because
which claims rest on which reading is itself part of the record.

### Verification pass against the Zotero PDFs (2026-08-21, later same day)

PDFs were subsequently added to Zotero, allowing the agent-reported findings to be
checked against local copies. Results:

| Claim | Status |
|---|---|
| Lin & Han: "Without refitting the nuisance estimator" | **CONFIRMED** in the SSRN copy (Zotero `TKPBYGYN`, 30pp). "refit" occurs **exactly once** in the whole document; no simulation section. arXiv and SSRN agree on the decisive point — **gap 4 closed.** |
| DiCiccio & Efron: no "percentile method" | **CONFIRMED** — "percentile method/interval/bootstrap" = **0 occurrences** in the JSTOR PDF (`MEIT8R7K`). "necessary and sufficient" present as quoted. |
| Zhan: KS-distance mechanism, not interval comparison | **CONFIRMED** (`XPEQV45A`, 18pp): "instruments are deemed weak iff the Kolmogorov-Smirnov (KS) distance below exceeds..." |
| Zhan: anti-pretest footnote removed in publication | **CONFIRMED** — `pre-test`, `pretest`, `Guggenberger` all return **0 hits** in the published version. |

At that point two gaps remained; **both have since been closed** — Zhan's appendix
arrived with the T&F supplementary package (which also carried the version-of-record
source, see §1), and the full DiCiccio & Efron discussion was added piece by piece
(§2). Nothing in this review now rests on an unread source.

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

### Verified against the version-of-record source + Online Appendix (2026-08-21)

The T&F supplementary package (`32399188.zip`) turned out to contain not only the
Online Appendix but `main.tex`/`main.pdf` — the version-of-record source. Every
finding above was re-checked against it, and all hold:

| Claim | Status against VoR / Appendix |
|---|---|
| KS-distance definition | **CONFIRMED**, `main.tex:194`: "deemed weak iff the Kolmogorov-Smirnov (KS) distance below exceeds 5\%" |
| Anti-pretest footnote not restored in the VoR | **CONFIRMED** — `pre-test`, `pretest`, `Guggenberger` all 0 hits in `main.tex`. The 2017 draft's caution genuinely did not survive into the published paper. |
| Theorem 4 is a pattern result, not consistency-of-detection | **CONFIRMED** by reading the proof (Appendix D): it derives how `mu^2` and `mu*^2` behave across the three regimes `0 <= delta < 1/2`, `delta = 1/2`, `delta > 1/2` (the last giving an asymptotic difference of `chi^2_k`). It is a concentration-parameter behavior result. There is **no theorem stating the test detects weak identification with probability tending to one.** |
| No ratios / DML in the appendix | **CONFIRMED** — `cross-fit`, `double/debiased`, `machine learning`, `semiparametric`: **0 hits**. The two `ratio`/`denominator` hits are incidental ("the concentration parameter divided by k"). |

**One refinement to the "no percentile interval" claim.** The paper *does* use Efron's
percentile interval — `main.tex:369`: "Such a confidence interval is called Efron's
percentile interval." But its **object is the KS statistic itself**, not the structural
parameter: the double bootstrap produces `KS**1,...,KS**B`, and their quantiles form a
CI *for KS*. There is still no percentile interval for `theta`, and still no
Wald-vs-percentile comparison of intervals for the estimand.

This makes the manuscript's misreading more understandable — "percentile interval" and
"weak identification" do co-occur in Zhan — but it does not rescue the citation. The
gauge manuscript attributes to Zhan a comparison between two intervals *for the
estimand*; Zhan builds one interval *for a distance measure*.

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

### The Hall & Martin discussion — and a correction to this review's synthesis

Obtained 2026-08-21 (pp. 212-214, JSTOR 2246111). Hall is the author of the Hall
(1988) Edgeworth comparisons that both DiCiccio & Efron and Owen (2025) build on, so
his view carries particular weight here. §5 (The Double Bootstrap):

> "One might summarize the respective theoretical drawbacks of percentile-t and BC_a
> methods by noting that the former are not transformation invariant, and **the latter
> are not monotone in coverage level**. As a utilitarian procedure **we favor a
> calibrated version of a simple method such as percentile. The percentile method is
> transformation respecting**; its calibrated form 'almost' respects transformations
> and is monotone in coverage level. **Also, it is not hindered by problems associated
> with ratios of random variables, which are sometimes the downfall of percentile-t.**"

**This partially contradicts the synthesis below, and the correction matters.** An
earlier version of this review concluded flatly that "the percentile bootstrap arm has
no citation support." That was too strong. In the published discussion of the very
paper the manuscript cites, Hall & Martin explicitly favor the percentile method and
name **ratios of random variables** as a case where it beats percentile-t.

Three scope conditions keep this from rescuing the manuscript's citation, but they
must be stated rather than used to wave the finding away:

1. **They favor a *calibrated* percentile**, not the raw quantile interval. Calibration
   is the double-bootstrap adjustment (Loh 1987; Beran 1987). `ward_residual()`
   implements raw percentile.
2. **The ratio remark contrasts percentile with *percentile-t***, not with BC_a, and
   the mechanism is percentile-t's dependence on a variance estimator — not a
   denominator approaching zero. It is not a statement about the near-null regime.
3. **The manuscript cites DiCiccio & Efron's main text**, which still does not present
   a percentile method at all. A discussant's endorsement is not the cited source.

**Two things this does change.**

- *For #32:* Hall supplies an additional, independent argument against BC_a —
  non-monotonicity in coverage level — from a source with no small-`n`-mean scope
  limitation, unlike Owen.
- *Constructively:* **calibration becomes a live third option** alongside the no-refit
  bootstrap and Fieller. The double bootstrap is normally prohibitive here (`B^2`
  full cross-fit refits), but under the no-refit reweighting scheme of Phase 2 it is
  `B^2` reweightings of an already-computed `phi` matrix. **Phase 2 is what makes
  Hall's recommended procedure affordable** — an argument for that phase that has
  nothing to do with citation hygiene.

### The rejoinder: DiCiccio & Efron pre-emptively rebut the manuscript's reading

Obtained 2026-08-21 (pp. 223-228). **This is the single most on-point passage in the
entire citation set**, and it is in the article the manuscript already cites (p. 226):

> "In fact it is difficult to run a good simulation study of confidence intervals
> methods. Besides the pitfalls mentioned earlier, and the cruel computational burden,
> there is the question of **interval length variability. One way to get better
> coverage accuracy is to make your intervals longer and more variable.** As an extreme
> example, one could choose `U` uniform on (0,1) and define [`theta[alpha] = 0` if
> `U < alpha`, `infinity` if `U > alpha`]. Then the interval `(-inf, theta[alpha])`
> would cover the true `theta` (or any other value) **with probability alpha**."

The authors are warning, in the same article, that **coverage attained by making
intervals longer is not evidence of anything** — and they give a degenerate construction
that hits nominal coverage exactly while carrying no information.

The gauge manuscript reports percentile-bootstrap coverage of 1.00 in 8 of 8 cells and
reads it as success: the bootstrap "restores coverage to nominal-or-above ... confirming
the ratio diagnosis." DiCiccio & Efron's rejoinder says in advance that this inference
is a trap. This is no longer a matter of a citation being used loosely — **the cited
source explicitly rebuts the conclusion drawn from it.**

It also independently vindicates A1's semantics. A1 flags intervals that are
uninformatively wide (flagged draws over-cover by +0.10, median 21x wider than
|truth|). That is precisely the failure mode D&E name, which is why A1 carries
information beyond width while A2 does not.

**On whether they accepted Hall's preference for calibrated percentile: partially.**
They endorse calibration as answering "how accurate are my confidence interval
coverages for my particular statistic and sample size?" — "The calibration methods of
Section 7 provide at least a partial answer" — but apply it to ABC/BCa rather than to
percentile. Agreement on the mechanism, not on the base method. They also concede BCa
and ABC "tend to be rather cautious improvements, sometimes not improving enough on
the standard intervals."

### The complete discussion (pp. 212-228): three positions, no winner

All remaining pieces obtained — Canty, Davison & Hinkley (214-219), Gleser (219-221),
Lee & Young (221-223). With Hall & Martin (212-214) and the rejoinder (223-228) the
discussion is now read in full, and the honest summary is that **the discussants do not
agree**:

| | Preferred method |
|---|---|
| Hall & Martin | calibrated **percentile** |
| Canty, Davison & Hinkley | **Studentized** (bootstrap-*t*), with a transformation |
| DiCiccio & Efron (rejoinder) | **BC_a / ABC**, plus calibration |

Hall & Martin say so explicitly in their opening — the paper "point[s] out that there
are no uniformly superior methods."

**The CDH position is a direct counterweight to the rejoinder passage quoted above, and
the two form an exchange worth reading together.** CDH (p. 218):

> "**Far from being a drawback**, in this problem the fact that the Studentized
> bootstrap method can give long confidence intervals **is precisely what gives it the
> best coverage** of the methods considered in our simulation study, and the
> 'conservativeness' of the BC_a method is what leads it to undercover."

DiCiccio & Efron then reply with the interval-length-variability warning and the
degenerate `(-inf, theta[alpha])` construction.

**Neither side is simply right, and the synthesis is the useful part.** CDH are correct
that length is not automatically a defect — it can be what legitimately buys coverage.
D&E are correct that length can also manufacture coverage vacuously. What follows is
not "long intervals are bad" but:

> **A coverage number is uninterpretable without an accompanying length or
> informativeness measure.**

That is the claim the gauge manuscript's Table 1 violates — it reports coverage alone,
in a design where the percentile arm hits 1.00 in every cell.

And our situation is diagnosable on exactly that axis rather than by assertion: the
flagged draws are median **21x wider than |truth|**, which places them in D&E's vacuous
regime, not CDH's earned one. A1 already separates the two.

One further CDH finding worth recording, since it cuts against BC_a independently of
Hall's non-monotonicity point (p. 217): in their study "percentile methods do very
poorly in the upper tail: the top endpoint of these intervals is too low. Unfortunately
**the same is true of the BC_a and the ABC methods, which do only as well as the much
simpler normal and basic bootstrap intervals.**"

Gleser's comment (219-221) concerns reproducibility — the "first law of applied
statistics," that a method should give the same answer on re-analysis, which Monte
Carlo methods violate — not ratios. Despite Gleser & Hwang (1987) appearing in the
rejoinder's bibliography, **no discussant raises the nonexistence of finite-length
confidence sets for ratio-type estimands**, and Fieller is never mentioned by anyone.
The connection to our near-null regime is ours to make, not something to cite them for.

### Verdict

**Premise supported; prescription inverted — but see the Hall & Martin correction
above.** DiCiccio & Efron's main text argues that being tail-aware means *estimating*
`z_0` and `a`, not reading raw quantiles, and the manuscript cites them for the rung of
the ladder they treat as deficient. What the discussion adds is that percentile is not
therefore indefensible — in *calibrated* form it has Hall's explicit endorsement, and
specifically for ratios.

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

**On Tang & Westling as the "refit case" citation — read directly 2026-08-21
(arXiv:2404.03064v2, 142pp), and it does NOT rescue the refit procedure.**

It was initially recorded here as the correct citation for a refit-per-resample
bootstrap. That framing was too generous. Their condition (B2) is high-level, so a
refit *could* satisfy it — but the paper names our exact procedure as the case at
risk (§3.2):

> "If eta*_n is constructed in an **exactly analogous manner using the bootstrap
> data** as eta_n is constructed using the original data, the bootstrap data has
> **replicated observations**, and the method of constructing eta_n is **sensitive to
> ties** in the data, **(B2) may not be satisfied**. ... for this reason and others we
> **do not require** that eta*_n be constructed in an exactly analogous manner to
> eta_n, so these issues can be avoided. In particular, the simplest approach for
> constructing eta*_n is to define **eta*_n = eta_n**."

`ward_residual()` does precisely the flagged thing: `sample.int(n, n, replace = TRUE)`
(Efron multinomial resampling, hence heavy ties) followed by a full `.corner_fit()`
nuisance refit on the tied resample — "an exactly analogous manner using the bootstrap
data." Whether (B2) actually fails depends on the learner's tie-sensitivity, which is
an open empirical question for our nuisance fits, but this is the paper's named
failure mode, not a generic caveat.

Their recommended simplest construction is `eta*_n = eta_n` — **no refit**, the same
place Lin & Han's theorem lives.

Caveat on using them at all: for the empirical bootstrap their condition leans on a
`P_0`-Donsker requirement (Gine & Zinn 1990), which cross-fitting exists to avoid;
§3.3 gives further sufficient conditions for smooth bootstrap sampling distributions.

**Net: both available theory papers point away from refit-per-resample and toward
holding the cross-fitted nuisances fixed.** There is no citation in this set that
endorses what the code currently does.

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

**The percentile bootstrap arm, as implemented, has no citation support in this
manuscript.** Its three supporting references either argue for something else
(DiCiccio & Efron's main text), cover a different procedure (Lin & Han), or address a
different question in a different model class (Zhan).

**Qualification, after reading the full discussion (see §2).** This is a claim about
*raw* percentile and about *these citations* — not a claim that the percentile method
is indefensible, and **not** a claim that some other method is established as correct
here. The discussion of the very paper the manuscript cites contains three
incompatible recommendations (Hall & Martin: calibrated percentile; Canty, Davison &
Hinkley: Studentized; DiCiccio & Efron: BC_a/ABC with calibration), and states outright
that "there are no uniformly superior methods."

So the defensible conclusion is narrower and more useful than "percentile is wrong":

1. **The manuscript asserts a resolution the literature does not support.** Its prose
   presents the percentile bootstrap as *the* appropriate construction for a skewed
   ratio, citing a paper whose own discussion could not agree on a winner.
2. **It reads a coverage number without a length check** — the specific error D&E's
   rejoinder warns against, and the one thing all sides of that exchange implicitly
   agree matters.
3. **Its bootstrap is separately unlicensed** by the DML citation (Lin & Han), which is
   a procedural defect independent of which interval type is best.

Points 2 and 3 stand regardless of how the method debate is settled. That is what makes
them actionable.

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

1. ~~**DiCiccio & Efron's discussion**~~ — **CLOSED.** All five pieces read (Hall & Martin 212-214; Canty, Davison & Hinkley 214-219; Gleser 219-221; Lee & Young 221-223; Rejoinder 223-228). See §2.
   (pp. 212-214, JSTOR 2246111) obtained and read**; it materially qualified this
   review's synthesis (see §2). Still unread: the remaining comments (pp. 214-228) and
   DiCiccio & Efron's rejoinder. Lower priority now that the most authoritative
   discussant is in hand, but the rejoinder would show whether the authors accepted
   Hall's preference for calibrated percentile.
2. ~~**Zhan's Supplementary Appendix A-F**~~ — **CLOSED.** Obtained via the T&F
   supplementary package, which also carried the version-of-record source
   (`main.tex`). All four claims re-verified; see the table in §1. Theorem 4's proof
   confirms there is no consistency-of-detection theorem, and the appendix contains no
   ratio, Fieller, or DML content. Extracted at scratchpad `zhan-suppl/`; the zip is
   at `~/Downloads/32399188.zip`.
3. ~~**Tang & Westling** (arXiv:2404.03064)~~ — **CLOSED.** Downloaded and read
   directly (v2, 142pp). Outcome reversed the plan: it does not license the refit
   procedure but names it as the at-risk case, and recommends `eta*_n = eta_n`. See §3
   above. PDF at `~/Downloads/Tang_Westling_2024_bootstrap_consistency_ML.pdf`, not
   yet in Zotero.
4. ~~Whether the arXiv and SSRN versions of Lin & Han are identical~~ — **CLOSED.**
   The decisive no-refit sentence, the single occurrence of "refit", and the absence
   of a simulation section all reproduce in the SSRN copy.
