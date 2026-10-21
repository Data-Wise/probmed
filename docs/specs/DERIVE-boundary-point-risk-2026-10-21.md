# Boundary point-estimation risk: derivation (T1) and check (T2)

Plan: `tasks/plan-boundary-point-risk.md`. Research only. No gate, estimator, threshold,
PR status or production code changes. Evidence: archived study on PR #55
(`docs/reports/VALIDATE-wasserstein-boundary-selection-2026-10-20.md`,
`data/wasserstein-selection-2026-10-20/replicates.csv`). No fresh datasets were generated.

## Setup

Binary cells. Contrasts `d = (d1, d2)` (indirect, direct), estimate `d_hat ~ N(delta, V)`
approximately, `s_j = sqrt(V_jj)`, `mu_j = delta_j / s_j`. Target
`P(d) = sqrt|d1| / (sqrt|d1| + sqrt|d2|)`. Estimators: raw plug-in, `k4` (plug-in with a
damped curvature shift), and hard selection: if exactly one `|d_hat_j| / s_hat_j <= lambda`,
snap to the axis (`P = 0` if j = 1, `P = 1` if j = 2).

## Derivation

**D1. Cusp bias is of order n^(-1/4), not n^(-1/2).** `sqrt|.|` is not differentiable at 0, so
the delta method does not apply near an axis. For `X ~ N(0, s^2)`,
`E sqrt|X| = s^(1/2) * E|Z|^(1/2)` with `E|Z|^(1/2) = 2^(1/4) Gamma(3/4) / sqrt(pi) ~ 0.822`.
When `delta_j = 0` the plug-in is therefore biased at rate `s_j^(1/2)`, and the other axis
enters through `sqrt|d_k|` at its true value. Ratio-of-means heuristic (relabeled after the independent audit; it replaces `E sqrt|X2|` by its mean inside the ratio and is not a first-order expansion):

- zero direct (cell 5): `d1 = .1123`, `s2 = .0428`: `P_hat ~ .3351 / (.3351 + .822*sqrt(.0428)) = .663`
  against truth 1, bias about -.337; archived raw bias -.326. The literal first-order term
  `-E|Z|^(1/2) sqrt(s2) / sqrt(d1)` is -.508 here and the quadrature bias at this s is -.325, so
  the heuristic is close partly by luck (the first-order term is only 64% of the bias at
  s = .0428 and reaches 96% at s/256).
- zero indirect (cell 6): `d2 = .1191`, `s1 = .01726`: `P_hat ~ .108 / (.108 + .345) = .238`
  against truth 0; archived raw bias +.2377 (literal first-order term +.313).

**D2. Snapping decomposes the bias.** Write `S` for the snap event and `c in {0, 1}` for the
axis value. Then `P_tilde = k4 + (c - k4) S` and
`bias = bias(k4) + E[(c - k4) S]`. Under approximate independence
`Pr(snap axis 2) = Pr(|Z + mu2| <= lambda) * Pr(|Z + mu1| > lambda)`, with
`Pr(|Z + mu| <= lambda) = Phi(lambda - mu) - Phi(-lambda - mu)`. Cell 2:
`mu1 = 5.43`, `mu2 = 2.73`, so `Pr(snap) ~ Phi(-.77) = .22`; archived `.2295`. Snapping
sends an estimate to 1 when the truth is .51, hence the +.096 bias.

**D3. Trade-off: no threshold separates "zero" from "small".** Raising lambda lowers the
activation probability at an exact zero (`Pr(|Z| > lambda)`; .05 to .005 for lambda 1.96 to
2.807) but raises `Pr(snap)` for small nonzero signals, where the snap error is
`|c - P(delta)|`. Near `delta_2 = c * s_2` with `d1` fixed,
`1 - P = sqrt(c s_2) / (sqrt(d1) + sqrt(c s_2))`, which tends to 0 only like `s_2^(1/2)`.
So the worst-case snap bias of a hard threshold also decays like n^(-1/4): thresholding
moves the n^(-1/4) bias from the exact null to a neighborhood of it; it does not improve
the rate. Weak-direct cell 7: `delta_2 = .0098`, `mu_2 = .23`, `1 - P = .228`, and any
tested threshold snaps to 1 with high probability (formula: `Phi(lambda - .23) - Phi(-lambda - .23)`
= .944 at lambda 1.96 and .994 at 2.807; archived .95 and .997), giving bias about +.20 to
+.23 (archived +.203, +.226). The probability is lower for smaller thresholds (.67 at
lambda = 1), so the claim holds for lambda of about 1.96 and above, not for "any lambda >= 1"
(corrected after Codex Review A, finding 1). Stronger selection helped cell 5
(bias/SD .243 to .076) and hurt cell 7 (1.87 to 7.88) for the same reason.

**D3'. Two-point lemma (promotes D3 for the local-minimax rate; derived by the literature
subagent and re-checked by me, not yet independently audited).** Model: `X ~ N(delta, s^2)`,
`d1 > 0` known, target `P(delta) = sqrt(d1) / (sqrt(d1) + sqrt|delta|)`. Take `delta_0 = 0`
and `delta_1 = s`. Then `P(delta_0) = 1` and `1 - P(delta_1) = sqrt(s) / (sqrt(d1) + sqrt(s))`,
which is of order `s^(1/2) / sqrt(d1)`. `TV(N(0, s^2), N(s, s^2)) = 2 Phi(1/2) - 1 = .383`, free
of `n`. Le Cam's two-point bound gives, for every estimator (any threshold or `lambda_n`
included), `max over {delta_0, delta_1} of E|T - P| >= (gap / 2)(1 - TV)`, so the local
minimax risk is at least a constant times `n^(-1/4)`. **Scope:** this is a local-minimax
statement over a `1/sqrt(n)` neighborhood. Pointwise at exactly `delta = 0` a constant (or a
threshold) has zero risk (Hodges), which is why the Hodges-problem decision (Q2) needs a
ladder. It says nothing about RMSE scaling at regular cells, where the rate is `n^(-1/2)`
(independent audit, Q5), so it does not license an `n^(-1/4)` tolerance there.

**D4. Interval and point disagree under strong selection.** The Fieller interval is
computed from the unselected `(d_hat, V)`; a snapped point can lie outside it
(archived point-in-interval .70 to .97 at lambda 2.807). This is an observation, not a
derived result: the snapped value is 0 or 1 regardless of where the inversion falls.

**Not established here.** (i) A formal minimax lower bound for all thresholding rules
(D3 is a local-rate argument, not a theorem). (ii) Whether any estimator, smooth or not,
passes all frozen gates in all seven cells. (iii) Non-existence of an unbiased estimator of
`sqrt|delta|` (suspected from non-analyticity; must be checked in the literature in T3
before it is used).

## T2 check against the archived replicates

Script: `scripts/check-boundary-risk-derivation.R`; output:
`docs/specs/data/check-boundary-risk-derivation-2026-10-21.txt`. Prediction: Gaussian
experiment `d_hat ~ N(true d, mean archived V)` per cell, 4e5 antithetic draws, arbitrary
seed 91 (not a study seed). The reimplemented point rules reproduce 900 archived rows to
1.2e-15. Gates stated before the run: bias within max(4 MCSE, .03); predicted/observed SD
in [.85, 1.15]; snap rate within .03.

| Gate | Result |
|---|---|
| Bias, 20 cell-method pairs | 20/20 pass |
| Snap rate, 10 pairs | 10/10 pass |
| SD ratio, 20 pairs | **17/20 pass** |
| Control A (delta-method bias = 0) | rejected in cells 5, 6, 7 as required |
| Control B (lambda - 1 at cell 2) | predicted .0387 vs observed .2295, rejected |

**Overall: FAIL by the pre-stated rule.** The SD misses are cell 4 select95 (predicted
.0905 vs observed .0770), and cells 5 and 7 select995 (.0620 vs .0378 and .0671 vs .0287).
The gate was not relaxed.

Interpretation (hypothesis, not tested): the failures are in near-saturated snap cells,
where SD is governed by the rare non-snap rate (predicted 1.7% vs observed .6% in cell 5;
SD scales with its square root). The rule uses studentized `|d_hat_j| / s_hat_j`, and in
the archive `cor(d_hat_1, v_hat_11)` is .72 to .77, so the lower tail of the studentized
statistic is lighter than the fixed-`V` Gaussian experiment assumes. A prediction that
resamples archived `(d_hat, v_hat)` pairs would test this; it is not done.

## T2 follow-up: studentization test (2026-10-21)

Script: `scripts/check-boundary-risk-studentization.R`; output:
`docs/specs/data/check-boundary-risk-studentization-2026-10-21.txt`. Same SD gate
([.85, 1.15]); 2e5 resamples; arbitrary seed 92; no fresh datasets.

- **Arm P (requested): resample archived (d_hat, v_hat) pairs.** SD gate 10/10, including the
  3 original misses. This arm cannot fail: each archived estimate is a deterministic
  function of its own pair, so a pair bootstrap reproduces the observed SD. It is not
  evidence for the mechanism.
- **Arm S (control): same rows, v_hat taken from a different row** (breaks the
  d_hat/v_hat dependence, keeps both marginals). SD gate 8/10. It **reproduces** the
  misses in cells 5 and 7 select995 (ratios 1.51 and 1.92; the fixed-V Gaussian gave
  1.64 and 2.34) and its snap rates (.984, .986) match the Gaussian prediction (.983,
  .982) rather than the observed (.994, .997). It does **not** reproduce the cell 4
  select95 miss (ratio 1.04 vs the Gaussian's 1.175).

Reading: the d_hat/v_hat dependence explains the two select995 misses (cells 5, 7), where
cor(d_hat_1, v_hat_11) is .72 to .78 and almost every fit snaps. It does not explain the
cell 4 select95 miss, whose cause is not identified here (the shuffle keeps the
marginals, so a Gaussian-shape or fixed-V assumption is a candidate, untested). This is
an attribution from one control, not a proof; both arms use the archived estimates.
The original T2 verdict (FAIL, 17/20) stands as recorded. Practical consequence for T3:
any SE or risk prediction for a snapped or near-snapped estimator must carry the
joint (d_hat, v_hat) dependence, not a fixed `V`.

## T2 follow-up 2: the cell 4 select95 miss (2026-10-21)

Script: `scripts/check-boundary-risk-cell4.R`; output:
`docs/specs/data/check-boundary-risk-cell4-2026-10-21.txt`. Same SD gate. Observed SD .0770
from 2000 replicates, 1.9% snapped (about 38 fits).

| Variant (one ingredient moved from the Gaussian toward the archive) | Predicted SD | Ratio | In gate |
|---|---:|---:|---|
| G0 Gaussian d_hat, fixed Vbar (the original prediction) | .0901 | 1.17 | no |
| G1 archived d_hat, fixed Vbar | .0788 | 1.02 | yes |
| G2 Gaussian d_hat, archived v_hat | .0904 | 1.17 | no |
| G3 Gaussian d_hat with the empirical Cov(d_hat), fixed Vbar | .0867 | 1.13 | yes (barely) |
| G4 archived d_hat, independent archived v_hat | .0797 | 1.04 | yes |

- **v_hat variability is not the cause:** G2 equals G0.
- **The d_hat distribution is:** archived d_hat has sd 1.010 / .978 times sqrt(Vbar diag),
  d1 is right-skewed (.25, excess kurtosis .17), and d2's mean sits +.0019 above truth.
  The snap event lies in the far tail (`z2 = .1191 / .0309 = 3.85` standard errors from the
  axis), where a 2% scale change or mild skew moves the tail probability a lot (snap rate
  .0286 predicted vs .0190 observed; G3 alone .0261).
- **Not sampling noise:** under G0, 4000 studies of n = 2000 give a study SD of
  .0808 / .0901 / .0994 (2.5 / 50 / 97.5%) and `P(SD <= .0770) = .0032`; 99.5% of studies fall
  inside the gate. So the Gaussian experiment is misspecified in this tail, not merely noisy.
- **Limits:** G1 and G4 resample the same 2000 d_hat that produced the observed SD, so they
  partly assume what they test; G3 resolves only part of the gap. The split between scale,
  skew and the mean offset is not isolated.

Verdict: the miss is best explained as tail sensitivity of a rare-snap SD to the
non-Gaussian, slightly narrower d_hat. That is a **supported hypothesis, not an identified
mechanism** (Codex Review A, finding 3): the ablations do not separate scale, skew and
mean offset, and G1/G4 reuse the archived d_hat. The original T2 FAIL (17/20) stands. For
the other two misses, the d_hat/v_hat dependence (cells 5, 7 select995) is the better
supported explanation. Bias and snap rates match, but note below what that does and does
not test.

**What T2 tests (Codex Review A, finding 2).** The scripts apply a copied point rule to
Gaussian or archived draws. Matching bias and snap rates therefore checks that the rule
is reproduced and that the Gaussian approximation is adequate for them. It does **not**
test D1's n^(-1/4) rate across sample sizes, D2's independence-product snap formula
(only its numerical value in cell 2 was compared by hand), or a squared-risk
decomposition. Those remain open (task T2b).

## T2b: testing the derivation itself (2026-10-21)

Script `scripts/check-boundary-risk-rate.R`; output
`docs/specs/data/check-boundary-risk-rate-2026-10-21.txt`. Tolerances are in the script
header and were fixed before the first run; none was changed. The only edit after the first
run removed a dead line (identical output).

| Gate | Result |
|---|---|
| A1: log-log slope of the raw bias at an exact zero contrast, 4 smallest s, in [.45, .55] | **.470, pass** (the first run used the wrong index and gave .459; corrected after the independent audit) |
| A2: bias / leading term `-E|Z|^(1/2) sqrt(s2) / sqrt(d1)` at the smallest s in [.93, 1.07] | **.964, pass** |
| A3 control: a regular configuration must NOT show slope .5 | slope 1.913, pass (the check discriminates) |
| B1: independence-product snap probability within .01 of the exact correlated value | 10/10 pass |
| B2: exact correlated value within 4 MC SE of a 4e6-draw simulation | 10/10 pass |
| B3: exact value within .03 of the archived snap rate | 10/10 pass |
| C1: branch-decomposition RMSE within [.85, 1.15] of observed | **8/10, FAIL** |

**Overall T2b: FAIL on C1 only.**

- **Rate (D1) is supported.** Local slopes over successive halvings are .445, .460, .471,
  .479, rising toward .5 as O(s) corrections shrink; the leading coefficient matches to 4%.
  The smooth (regular) case has slope about 2, so the bias at a zero contrast is genuinely
  of order n^(-1/4) in this model (independent Gaussian contrasts, deterministic quadrature).
  The slope gate passed narrowly because the asymptote is approached from below; the gate
  was not adjusted.
- **Snap formula (D2) is supported.** The independence product agrees with the exact
  correlated probability to better than .004 in every cell here (correlation is about
  -.33), so the independence simplification is harmless at these cells; the exact value
  agrees with simulation and with the archived rates.
- **Risk decomposition (C1) fails in the two rare-snap cells already flagged:** cell 4
  select95 (predicted RMSE .0910 vs observed .0772) and cell 5 select995 (.0621 vs .0379).
  These are the cells where the Gaussian SD predictions also missed, so C1 does not
  independently confirm the tail-shape hypothesis; it shows the same limitation. Cell 7
  select995, an earlier SD miss, now passes because its RMSE is bias-dominated.
- **What this does not show:** T2b tests D1 and D2, not D3. D3's rate statement is now backed by the two-point lemma D3' (local minimax only, pending independent audit).

## Consequence for T3 (if the author approves continuing)

The derivation implies a candidate must (a) correct the n^(-1/4) cusp bias without a
discontinuity, (b) keep regular and weak-direct cells at least as good as `k4`, and (c) get
its SE from something that respects studentization. D1 and D3 do not say such an estimator
exists; a negative result remains a valid outcome.
