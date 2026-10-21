# T3: gate feasibility and a smooth candidate (redrafted 2026-10-21)

Plan: `tasks/plan-boundary-point-risk.md`. Design only: no fresh datasets, no estimator frozen,
no gate changed. Inherits Checkpoint A with three limitations: (1) T2b C1 fails 8/10 in
rare-snap cells; (2) D3 (no threshold improves the n^(-1/4) rate) is an argument, not a theorem;
(3) the cell 4 mechanism is a hypothesis. Reviewed by Codex (B1); its three findings are folded in
below and triaged in `REVIEW-codex-boundary-point-risk-B1-2026-10-21.md`.

## 1. Result

At n = 500, **no [0,1]-valued estimator can pass the frozen `|bias|/SD < .10` gate in both cell 5
and cell 7 at the population level**: the two data laws are too close (TV .140) for any estimator
to separate their means by the .179 the gate requires. For the frozen empirical gate (2000
replicates per cell) a pass is possible only by sampling accident: probability <= 8.3e-5
(Bernstein, worst case over the population means; 7e-7 by normal approximation). A smooth candidate that
reduces stress-cell bias is specified below; it fails the complete gate in every cell it targets.
T3 therefore does not recommend freezing any candidate against the gate as written.

## 2. The feasibility argument

Cells 5 (direct = 0, P = 1) and 7 (direct = .05, P = .7721) share alpha, interaction and
covariates; their data laws differ only in the `Y | A, M, C` logit by `.05 * A`. For any
estimator `T` with values in [0,1] (randomized allowed), writing `E_k` for its mean in cell k:

1. Cell 5: `bias = 1 - E5`, and `SD5 <= sqrt(E5 (1 - E5))`. The gate `bias/SD < .10` then
   requires `E5 >= .9901`.
2. Cell 7: `|E7 - .7721| < .10 SD7` with `SD7 <= sqrt(E7 (1 - E7))` requires `E7 <= .811`
   (a looser bound, using `SD <= .5`, gives `.8221`).
3. `|E5 - E7| <= TV(P5, P7)`, so passing both needs `TV >= .9901 - .811 = .179`
   (`.168` under the loose bound).
4. **TV estimate** (`scripts/check-gate-feasibility-tv.R`, exact Bernoulli likelihood ratio,
   20,000 datasets per shift): **.1404 (MC SE .0012)** for the actual .05 shift. Controls:
   shift 0 gives 0, shift .25 gives .623, shift .5 gives .921. Output:
   `docs/specs/data/check-gate-feasibility-tv-2026-10-21.txt`.

**Empirical gate.** For values in [0,1] the sample SD obeys `s^2 <= R/(R-1) m (1 - m)`
deterministically, so steps 1-3 hold for the sample means (approximately, since `R/(R-1) =
1.0005`): an empirical pass needs `mbar5 >= .9901` AND `mbar7 <= .8112` separately. The two
cells use independent datasets, so `P(pass) <= P(mbar5 >= .9901) * P(mbar7 <= .8112)`, with
population means constrained by `|E5 - E7| <= .140`. Maximizing over `E5` (taking `E7 = E5 -
TV`, the most favorable value), Bernstein (variance `E(1-E)/R`, range `1/R`) gives
**<= 8.3e-5** and the normal approximation **<= 7e-7**
(`scripts/check-gate-feasibility-tv.R`). An earlier version of this paragraph bounded a
difference of means with a non-worst-case variance (.0114) and quoted .0042; the
independent audit showed that derivation was wrong (a worst-case difference variance gives
about .007), although the conclusion held. It is replaced by the product bound above.

**Assumptions:** estimator values in [0,1]; the gate applies to both cells as frozen; the two
cells' replicates are independent; TV is a Monte Carlo estimate, not a closed form.

**Why the [0,1] assumption carries the result (literature subagent, 2026-10-21; its own
derivation, not yet independently checked).** The two-point inequality gives
`(E1 T - E0 T)^2 <= chi2 * Var0(T)`, so it only bounds the mean gap relative to the SD.
Bias/SD is scale-free: an estimator allowed outside [0,1] can pass the gate by inflating its
variance, which is the Doss-Sethuraman / Hirano-Porter mechanism (small bias forces large
variance; Hirano and Porter 2012 abstract read by the subagent: "no locally asymptotically
unbiased estimators and no regular estimators"). My bound `SD <= sqrt(E(1-E))` is what blocks
this inside [0,1]. Consequence: bias/SD alone is a gameable gate, which supports gating RMSE or
risk, not bias/SD.

**Reading.** This does not say the Fieller interval work is wrong (its coverage passes), nor
that the gate was badly meant. Truths .228 apart differ by only about .23 standard errors in
`d_2`, and bias/SD cannot be small at both points.

## 3. A smooth candidate `C(kappa, tau)` (binary outcome, no threshold)

`psi_j = ((d_j^2 - kappa v_jj) + sqrt((d_j^2 - kappa v_jj)^2 + (tau v_jj)^2)) / 2`,
`m_j = psi_j^(1/4)`, `P_hat = m_1 / (m_1 + m_2)`.

- Smooth and strictly positive for `tau > 0`, defined at `d = 0`; `psi` is homogeneous of
  degree 2 in `(d, s)`, so `P_hat` is scale-invariant; swapping `d_1`, `d_2` gives `1 - P_hat`.
- `kappa = 1` makes `d^2 - v` unbiased for `delta^2` when `v` is fixed or independent of `d`; in the archive `cor(d_hat_1, v_hat_11)` is .72 to .78, so this is not exact there. Larger `kappa` shrinks harder toward 0.
- It keeps the n^(-1/4) rate; only the coefficient `c(kappa, tau)` replaces .822 (D1, D3).

**Predicted bias** (Gaussian experiment at the archived cell parameters; bias was validated
20/20 in T2; `scripts/explore-smooth-boundary-grid.R`, output
`docs/specs/data/explore-smooth-boundary-grid-2026-10-21.txt`; cells 2, 4 regular, 5 zero
direct, 6 zero indirect, 7 weak direct):

| kappa, tau | c | cell 2 | cell 4 | cell 5 | cell 6 | cell 7 |
|---|---:|---:|---:|---:|---:|---:|
| 0, .5 (near plug-in) | .94 | +.007 | +.003 | -.365 | +.268 | -.139 |
| 1, .5 | .70 | +.020 | +.008 | -.298 | +.223 | -.073 |
| 2, .5 | .56 | +.037 | +.014 | -.255 | +.198 | -.030 |
| 3.84, .5 | .42 | +.074 | +.028 | -.213 | +.186 | +.013 |
| k4 (archived) | | +.004 | -.001 | -.333 | +.239 | -.108 |

Against k4, cell by cell, at `kappa = 3.84, tau = .5`: better absolute bias in cells 5, 6 and 7
(cell 7: +.013 vs -.108, a sign crossing from over-shrinkage, not accuracy: bias/SD is still
.22), worse in the regular cells 2 and 4 (+.074 and +.028 vs +.004 and -.001). Within the
family (against `kappa = 0, tau = .5`): the zero-contrast coefficient falls 2.2x (.94 to .42;
1.96x against the raw .822), the cell 5 bias falls 1.71x and the cell 6 bias 1.44x, while
the regular-cell bias grows about 11x (cell 2) and 8x (cell 4). Against k4 the regular-cell
growth is about 18x and 28x.

**Pass/fail expectation (T3 verifiability criterion).** Predicted `|bias|/SD` is 2 to 11 in
cells 5 and 6 at every grid point, and .22 in cell 7 even at the best grid point (gate .10);
the regular cells fail for `kappa >= 1` (.17 to .58). So the family improves stress-cell bias
and still fails the complete gate where it was designed to help. SD here is from the Gaussian
experiment and is indicative only (T2: unreliable in rare-snap cells).

## 4. Literature (what was and was not checked)

- **Derived here, not cited:** for `X ~ N(delta, s^2)`, `E_delta g(X)` is real-analytic in
  `delta` for integrable `g`, while `sqrt|delta|` is not analytic at 0, so no unbiased
  estimator of `sqrt|delta|` exists (same for `|delta|`).
- **Doss and Sethuraman (1989), Ann. Statist. 17(1), 440-442:** a search summary says bias
  reduction toward zero forces variance up when no unbiased estimator exists. The paper could
  not be opened, so this is a pointer only.
- No source was found for plug-in bias-reduction rates. Nothing here depends on a cited result.

## 5. Options for the author (no gate is changed here)

1. **Do not freeze (T5) or run (T6) a candidate against the unchanged all-cells bias gate.**
2. **Redesign the gates and pre-register them before any fresh data** (the author's proposal is
   pending). Check each proposed gate against the same two-point condition: for nearby cells,
   feasibility needs TV between their data laws to exceed the truth gap the gates would force
   apart. The change must be justified independently of which candidate it favors and cannot
   be applied retroactively to archived results.
3. **Record the negative-result endpoint:** a derived rate (D1), a derived trade-off
   (D2/D3), a measured infeasibility of the frozen gate (section 2), and the Fieller interval
   as the supported deliverable.

T4 and T5 are held. Review B stays scheduled before any freeze.
