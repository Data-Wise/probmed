# Correction sweep: propagating the 2026-08-22 refutation to every surface

**Created:** 2026-08-22
**Trigger:** `REVIEW-2026-08-22-adversarial-refutation.md`
**Purpose:** the refuted claims were written to six surfaces over one session. This
inventories them and specifies the correction for each, so the retraction is neither
incomplete nor overreaching.

---

## The governing distinction — read before editing anything

The review measured **2 of 8 cells**, both `tint = 0`, continuous-Y, at
`oe_snr` median **8.2 / 14.9**. Those are the *strong-identification* cells. It did not
test binary-Y, `tint = 0.8`, or anything near the null.

Three categories, and conflating them would replace one error with another:

| | Claim | Status |
|---|---|---|
| **REFUTED** | Mechanism claims (heavy tails cause over-coverage; reweighting `phi` is a fix; `oe_snr = z` partitions the space) | Wrong as reasoning, at any `n`. Retract outright. |
| **REFUTED IN THE TESTED REGIME** | "The percentile interval is vacuously wide" | False where measured (0.80-0.83x a calibrated normal interval). **Untested near the null.** Scope the retraction; do not generalize it. |
| **UNAFFECTED** | What Zhan / DiCiccio & Efron / Lin & Han / Owen / Hall actually say | Source readings stand. Only their *application to our case* is in question. |

**The `21x wider than |truth|` figure is in the third-and-a-half category**: it comes
from A1-*flagged* draws in the 48k weak-ID grid, a regime the review never entered. It
is **not refuted** — it is unreplicated in the tested cells, which were strong-ID and
where intervals measured *narrower* than calibrated. Both can hold simultaneously. Do
not delete it; scope it.

---

## Surface-by-surface

### 1. Issue #11 comment — PARTIAL, urgent

The comment asserts two things. They now diverge:

- **A1 = Fieller boundedness criterion** — **SURVIVES** (C2). Keep.
- **"The bootstrap and the Fieller set partition the parameter space at a threshold A1
  already computes"** — **REFUTED** (C4). A3.2(iii) constrains the *population*
  Jacobian `J_0 = -OE`, a fixed constant satisfiable at every `n`; `oe_snr` is a
  *sample* statistic growing like `sqrt(n)`. Fixed-vs-local asymptotics conflation.

Also add the two weakenings the review found on the surviving half:
- the identity is contingent on `ci_level` (`fa` uses `zc`; A1 uses a fixed
  `oe_snr_threshold = 2` — at `ci_level = 0.99` Fieller needs `snr > 2.576`);
- under `reps > 1`, lines 302/366 use unadjusted `var(pOE)/n` while 275-276 add the
  fold-split component to `seW`/`seP` only, so **the Fieller set is declared bounded
  more often than the rest of the machinery justifies**;
- "A1 detects non-identification rather than narrowness" is *partly circular* — `|OE|`
  enters the interval width through the same `1/OE` factor.

**Priority: highest.** A public claim known to be wrong, same case that justified the
#32 correction.

### 2. Issue #32 — DONE (2026-08-22)

Retraction posted. Records that the diagnosis cuts *toward* BCa (shape is BCa's domain)
rather than away, and withdraws the "close as answered" recommendation. Nothing further
unless the `phi` decomposition changes the picture.

### 3. `.STATUS` item 4(a) — SCOPE, don't delete

The *observation* survives: the manuscript reports 1.00 coverage and frames it as
success without a length check. The *interpretation* does not:

- delete "1.00 in every cell ... is uninformative, not calibrated" — measured 0.80-0.83x
  calibrated width in the tested cells says otherwise;
- delete the over-dispersion speculation ("more consistent with an over-dispersed
  resample distribution");
- correct "1.00 in 8/8" to **6/8 at 1.000 and 2/8 at 0.9995** (`mcse = 5e-04` inverts to
  one miss each);
- keep the core finding — a coverage number reported without a length or
  informativeness measure is uninterpretable — since that is what both sides of the
  CDH/D&E exchange agree on and it does not depend on the refuted mechanism.

### 4. `SPEC-2026-08-21` — header done, body stale

The superseded header is in place. The body still asserts refuted content at:
`:65` (heavy tails as mechanism), `:80` ("merely widens until it vacuously covers"),
`:152-165` (second mechanism + the no-refit fix), `:211 / :330 / :406` (the partition),
`:395` (the D&E "longer and more variable" framing applied to our case).

**Do not patch it inline.** Rewrite it after the `phi` decomposition produces
measurements, so the replacement is written against numbers rather than another guess.
Until then the header does its job: the document is a record of what was proposed, and
is marked not-to-execute.

### 5. `LITREVIEW-2026-08-21` — separate sources from application

The literature findings are unaffected and should not be touched. What needs correcting
is where the review *applies* them to our case:

- `:248-250` — "calibration becomes a live third option ... Phase 2 is what makes it
  affordable." The enabling premise is refuted. Keep Hall's endorsement of calibrated
  percentile; drop the claim that our Phase 2 unlocks it.
- `:272` and `:320-330` — the rejoinder section applies "longer and more variable" and
  concludes our draws sit "in the vacuous regime, not CDH's earned one." **Scope this
  to the near-null flagged draws** and record that in the tested strong-ID cells the
  measurement runs the other way.
- `:525` — the synthesis' "covers at 1.00 in 8 of 8 cells." Correct to 6/8 + 2/8.
- `:544` — "the no-refit multiplier bootstrap is now doubly motivated." Refuted; strike.

### 6. `R/gauge-pmed.R` — two in-code corrections, both still owed

- **The Zhan attribution** (roxygen + the `weak_id` warning text). Unaffected by the
  refutation — A2 is not Zhan's statistic. Restate as an adaptation rather than
  dropping the citation outright (C6 weakened, not refuted).
- **New, from the review:** add a note that A2's ~2x null baseline is a consequence of
  `sigma-hat_an` being biased low by ~30% in all 8 cells. A diagnostic whose baseline
  is 2x *because its own denominator is biased* is partly measuring SE bias rather than
  identification strength. This subsumes the 97%-width-explained finding and is the
  better characterization.

### 7. The manuscript (cross-repo) — plan changed, nothing written yet

No corrections owed, because nothing was written. But the planned Phase 5 edits are now
partly wrong and must not be executed as specified: the "percentile reaches 1.00 by
being far wider than warranted" rewrite is refuted in the tested regime, and the D&E
rejoinder framing does not apply there.

**What survives for the manuscript, independent of all of this:**
- the **Lin & Han miscitation** (C3 survives) — the code refits, the theorem does not
  cover refitting;
- the **Zhan mechanism/setting misattribution** (C6 weakened but standing);
- the **DiCiccio & Efron clause re-scoping** — their main text presents no percentile
  method; the premise they support is that symmetric intervals mis-cover skewed
  estimands.

These three are citation-accuracy findings and do not depend on any coverage
measurement. They can proceed whenever the cross-repo go-ahead is given.

---

## Ordering

```
1. Issue #11 correction ............... ~10 min, unblocks nothing, stops a wrong claim
2. .STATUS item 4(a) scoping .......... ~15 min, in-repo
3. LITREVIEW application fixes ........ ~30 min, in-repo, leaves sources untouched
4. R/gauge-pmed.R comment corrections . ~20 min, needs a branch (code file)
   |
   +-- 5. Decompose the `phi` gap ..... the real work; produces the numbers
              |
              +-- 6. Rewrite SPEC against those numbers
                        |
                        +-- 7. Manuscript (cross-repo, needs go-ahead)
```

Steps 1-3 are documentation-only and can land on `dev` directly. Step 4 touches a code
file's comments, so it needs a branch under the local guard rules.

## The one measurement that governs everything downstream

**Decompose `seW_an_ratio = 0.65-0.83` into fold-split vs nuisance-refit vs remainder.**

It was visible in the spec's own Phase 0 table from the start and never pursued. Until
it is understood, every interval built from `phi` — Wald, any reweighted bootstrap, and
**Fieller** — is too narrow by the same factor, and the Fieller set is declared bounded
more often than it should be. `R/corner.R:29` (`folds <- sample(...)` per call) is one
identified source; the package already concedes the component is material by adding
`var(W_reps)/reps` at `reps > 1`.

**Extend it beyond the review's coverage before drawing conclusions**: the review ran
2 of 8 cells, both continuous-Y at `tint = 0`. The binary-Y cells carry the untested
`seW_bt_ratio` outliers (5.06 and 16.62), and the near-null regime is where A1, Fieller,
and the manuscript's actual claims all live.
