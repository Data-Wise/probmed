# Adversarial review of SPEC-2026-08-21: Phase 2 refuted, Phase 0 diagnosis wrong

**Date:** 2026-08-22
**Status:** the spec's primary fix is refuted by measurement; the spec must not be
executed as written
**Method:** independent fresh-context review with repo access and no sight of the
reasoning that produced the spec. The reviewer ran its own experiments rather than
arguing from the documents.

---

## Verdicts

| Claim | Verdict |
|---|---|
| C5 — no-refit bootstrap is the fix | **REFUTED** (measured) |
| C1 — heavy tails explain 1.00 coverage | **REFUTED** (mechanism backwards) |
| C4 — bootstrap/Fieller partition at `oe_snr = z` | **REFUTED** (category error) |
| C2 — A1 gate = Fieller boundedness | SURVIVES, two real weakenings |
| C6 — A2 is a Zhan misattribution | WEAKENED (worse defect found) |
| C3 — Lin & Han miscitation | SURVIVES (but see the fiat problem) |

## The experiment

Reproduced the estimator's internals (`probmed:::.corner_fit`, K=5), 40 reps x B=399,
cells `tint=0, binY=FALSE` at n=800 and n=3000 — the two **easiest** cells in the grid
(`oe_snr` median 8.2 / 14.9, nowhere near the near-null regime). Three arms on the same
fitted `phi`.

| within-rep ratio | n=800 | n=3000 |
|---|---|---|
| `sd(no-refit boot) / seW_an` | **1.046 ± 0.049** | **0.996** |
| `sd(refit boot) / seW_an` | 2.256 | 1.889 |
| refit resample kurtosis | 14.9 | 15.6 |
| refit `(q97.5-q2.5)/sd` | **3.858** | **3.756** (normal = 3.92) |
| coverage wald / refit / no-refit | .875 / 1.000 / **.800** | .825 / 1.000 / **.775** |
| mean width vs calibrated `2*1.96*empSD` | **0.80x** | **0.83x** |

---

## C5 — the proposed fix is the analytic SE in disguise

Reweighting a fixed `phi` reproduces the delta-method variance **identically**:
`th* = sum(phi_i w_i)/n` has variance `Sigma_phi/n`, which is exactly what
`R/gauge-pmed.R:270-272` computes. Measured at 1.046 / 0.996, every rep in [0.95, 1.17];
Gamma weights change nothing (1.041).

**The spec contradicts itself and never notices.** Phase 0 finding #1 states the
`phi`-based analytic SE under-disperses by 0.65-0.83x and concludes "nothing to change
here." Phase 2 then proposes replacing the bootstrap **with that same variance**.

Phase 0b's decision rule — "if coverage drops from 1.00 toward nominal, the refit was a
live contributor" — is therefore not a discriminating experiment. Any drop lands at the
analytic arm's level by construction.

Paired on the same 40 datasets: refit covers 40/40, no-refit 32/40 (n=800) and 31/40
(n=3000); all 8 discordant pairs go one way (exact p = 0.0078); against nominal 0.95,
p = 7e-4. **The fix overshoots and lands at or below the current Wald default, in the
two easiest cells in the grid.**

**Mechanism the spec missed, verified independently:** `R/corner.R:29` does
`folds <- sample(rep(1:K, length.out = n))` — **the fold split is redrawn on every
call.** So the refit bootstrap carries a fold-split component *plus* a nuisance-refit
component that fixed-`phi`, fixed-partition reweighting removes by construction. The
package already concedes this component is material: `R/gauge-pmed.R:274-277` adds
`var(W_reps)/reps` to `seW` when `reps > 1`. Phase 2's instruction to "keep the fold
partition fixed across draws, per §2.1" **deletes a variance source the point estimator
actually has.**

Better options, in order:
1. Keep the refit and fix the *reporting* — the coverage number may not be the problem
   (see C1).
2. If the refit must go, the honest replacement is no-refit reweighting **plus** an
   explicit fold-split/nuisance component (`reps > 1`, already implemented) — not
   `phi` alone.
3. `m`-out-of-`n` / subsampling, or a tie-robust learner, keeps the refit while
   addressing Tang & Westling's actual concern.

Calibrated percentile (the Phase 2 "unlock") calibrates *toward* nominal from a base arm
that would now sit at 0.78. Calibration repairs a shape error, not a missing variance
component.

## C1 — the diagnosis has the mechanism backwards

1. **Heavy tails compress the quantile range relative to SD** (t_3: 3.67 < 3.92).
   Measured: the refit resample distribution *is* heavy-tailed (kurtosis 15) — and its
   `(q97.5-q2.5)/sd` is **3.86 / 3.76 < 3.92**. The percentile interval is therefore
   *narrower* than +/-1.96 sd of the same replicates. Heavy tails work **against**
   over-coverage. "Correct SD + 1.00 coverage implies heavy tails" is a non-sequitur in
   the direction the spec needs.
2. **"Vacuously wide" is empirically false.** Mean refit width is **0.80x / 0.83x** a
   perfectly calibrated normal interval. It covers ~100% while being *narrower* than the
   interval you would draw knowing the true sampling SD. **Length is not what buys the
   coverage** — which removes the premise of the "#32 closes on the merits" argument and
   of the DiCiccio & Efron "intervals longer and more variable" framing in Phase 5.
3. **`seW_bt_ratio` is the wrong statistic** — a mean across reps of a within-rep SD,
   compared against the across-rep SD of the point estimate. The per-rep distribution is
   right-skewed (median/empSD 0.77 vs mean 0.83 at n=800). A ratio near 1 never licensed
   "each rep's bootstrap SD is correct."
4. **"Exactly 1.000 in 8/8" was a rounding artifact** — verified independently: the CSV
   stores `mcse_covW_pct = 5e-04` for two cells, which inverts to coverage **0.9995**,
   one genuine miss each. Six cells are 2000/2000, not eight. The spec argued from an
   exactness that is not in the data.

**What the reviewer measured instead:** the *sampling* distribution of `W-hat` is itself
heavy-tailed (median|err|/empSD = 0.39, 0.51 vs 0.674 normal), so `empSD` overstates the
typical error; and the percentile interval's **asymmetry** does the covering. At
n = 3000, `|W-hat - trW|` exceeds the symmetric half-width in **12.5%** of reps while
coverage is ~1.000 — the interval is **shifted, not wide**. The right diagnostic is the
distribution of `|W-hat - trW| / half-width`, not SD ratios or tail shape.

## C4 — the partition is a fixed-vs-sample category error

A3.2(iii) constrains the **population** Jacobian `J_0 = -OE`, a fixed constant. For any
DGP with `OE != 0` one can choose `c_0 < |OE| < c_1`, so the condition holds **at every
n**, independent of `oe_snr`. `oe_snr` is a **sample** statistic growing like `sqrt(n)`
that crosses `z` at some finite `n` for any `OE != 0`.

Identifying a fixed-parameter regularity condition with a sample-level test event
conflates fixed and local asymptotics. The defensible version would be a local-to-zero
argument (`OE = c/sqrt(n)`, so `oe_snr = O(1)`) — a weak-identification framing that
A3.2(iii) does not make. **"They partition the parameter space at the same threshold" is
numerology with a derivation pasted on.** The embedding `J_0 = -OE` is itself fine.

Boundary behavior finishes it: just above `z`, `fa -> 0+` sends the Fieller roots to
+/-infinity while the bootstrap is licensed only vacuously — the handoff switches
exactly where both arms degrade.

(The `c_1` upper-bound concern was a nothingburger: `OE` is a contrast of bounded corner
means, so `c_1` is non-binding.)

## C2 — the Fieller identity survives, with two weakenings

Arithmetically real: `seOE^2 = var(pOE)/n` (`:270`) and `VOE = var(pOE)/n` (`:366`) are
the same object, so `fa > 0 <=> |OE|/seOE > zc`.

- **It is a coincidence of two defaults, not structural.** `fa` uses `zc` derived from
  `ci_level`; A1 uses the fixed argument `oe_snr_threshold = 2`. At `ci_level = 0.99`
  Fieller needs `snr > 2.576` while A1 still fires at 2. Nothing couples them.
- **Under `reps > 1` the identity holds but coherence breaks.** Lines 302 and 366 both
  use unadjusted `var(pOE)/n`, while 275-276 add the fold-split component to
  `seW`/`seP` only. So Fieller and A1 ignore a variance component the Wald arm includes,
  and **the Fieller set is declared bounded more often than the rest of the machinery
  justifies.**

"A1 detects non-identification rather than narrowness" is *partly* circular: `|OE|`
enters the W-interval width through the same `1/OE` factor, so the +0.036 residual
measures what that particular width control failed to absorb, not a demonstrably
distinct mechanism.

## C3 — miscitation survives; the fiat around it does not

The quote is verified and the code demonstrably refits. But two things the spec must own:

- Its own Phase 2 records that Lin & Han average fold-wise *solutions* while
  `ward_residual` pools then divides — so Theorem 3.2 does not cover the no-refit
  version **either**. That contradicts Phase 5's "after Phase 2 the citation becomes
  correct as written."
- The spec resolves theory-vs-measurement by fiat: "correct on theoretical grounds
  **regardless of what the re-run measures**." The measurement now exists and says the
  refit carries real dispersion the fixed-`phi` construction cannot. **Choosing citation
  hygiene over measured calibration is the actual error.**

## C6 — attribution weakened; a better attack exists

A2 is a scale-only functional of the same bootstrap-vs-normal discrepancy Zhan's KS
measures — same family, cruder statistic, blind to shape and location. "Cite as an
adaptation, state it is not Zhan's statistic" is calibrated; "drop it" is too harsh.

The sharper defect: because `sigma-hat_an` is biased low by ~30% in **all 8 cells**
(0.65-0.83), the width ratio carries a large regime-independent offset — which is
exactly why the code needs a 3x threshold "calibrated above the ~2x baseline."
**A diagnostic whose null baseline is 2x because its own denominator is biased is
measuring SE bias, not identification strength.** That subsumes the 97%-width-explained
finding.

---

## The single most dangerous unexamined assumption

**That `phi` is a sufficient statistic for the sampling variability of `W`.**

Every downstream step assumes it: Phase 2 reweights `phi` alone; Phase 3's Fieller uses
`var(pOE)/n` from `phi` alone; Phase 4's positive control (`pct_oe_irreg`) is computed
from it; Phase 5 writes the manuscript on it.

The spec's own table already falsifies it — `seW_an_ratio = 0.65-0.83` in 8/8 cells —
and the measurement pins the gap at `sd_refit/sd_no-refit` = 1.9-2.3, with the
per-call fold redraw as one identified source.

**Until that gap is decomposed (fold-split vs nuisance-refit vs remainder) and repaired,
every interval built from `phi` — Wald, no-refit bootstrap, and Fieller — is too narrow
by the same factor, and the Fieller set will be declared bounded more often than it
should be.** That is the shared failure mode of Phases 2, 3, 4 and 5, and no phase in
the spec is scoped to find it.

---

## Caveats on this review, stated by the reviewer

- **nrep = 40 cannot distinguish coverage 1.000 from 0.96** (MCSE 0.047). The width and
  SD-ratio comparisons are paired within-rep and carry the argument; the coverage column
  is directional confirmation.
- **Only 2 of 8 cells were run**, both continuous-Y at `tint = 0`. They are the easiest
  cells, which is why undercoverage there is damning for C5 — but the binary-Y and
  near-null cells are untested.

## Independently verified before acceptance

- `R/corner.R:29` redraws folds per call — confirmed by direct read.
- `mcse = 5e-04` inverts to coverage 0.9995 at nsim = 2000 — confirmed by computation
  against the stored CSV values (`0,0,5e-04,5e-04,0,0,0,0`).

## What must happen before any code is written

1. **Do not execute Phase 2 as written.**
2. **Decompose the `seW_an_ratio = 0.65-0.83` gap** — fold-split vs nuisance-refit vs
   remainder. This is the real defect and it was hiding in plain sight in the spec's own
   Phase 0 table.
3. **Re-derive the C1 diagnosis** using `|W-hat - trW| / half-width` and interval
   asymmetry rather than SD ratios.
4. **Withdraw the "vacuously wide" framing** from Phase 5 and from the issue #32
   comments until the width measurement is confirmed at more cells — the intervals
   measured *narrower* than calibrated, not wider.
