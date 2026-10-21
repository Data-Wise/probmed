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
replicates per cell) a pass is possible only by sampling accident: probability <= .0042
(Bernstein, worst-case variance; .0003 by normal approximation). A smooth candidate that
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
1.0005`). With `E5 - E7 <= .140`, an empirical pass needs sampling noise of at least .039 in
`mbar5 - mbar7`. Its worst-case SD is .0114, so the needed noise is 3.4 SD, giving the bounds
above for any single estimator on independent datasets.

**Assumptions:** estimator values in [0,1]; the gate applies to both cells as frozen; the two
cells' replicates are independent; TV is a Monte Carlo estimate, not a closed form.

**Reading.** This does not say the Fieller interval work is wrong (its coverage passes), nor
that the gate was badly meant. Truths .228 apart differ by only about .23 standard errors in
`d_2`, and bias/SD cannot be small at both points.

## 3. A smooth candidate `C(kappa, tau)` (binary outcome, no threshold)

`psi_j = ((d_j^2 - kappa v_jj) + sqrt((d_j^2 - kappa v_jj)^2 + (tau v_jj)^2)) / 2`,
`m_j = psi_j^(1/4)`, `P_hat = m_1 / (m_1 + m_2)`.

- Smooth and strictly positive for `tau > 0`, defined at `d = 0`; `psi` is homogeneous of
  degree 2 in `(d, s)`, so `P_hat` is scale-invariant; swapping `d_1`, `d_2` gives `1 - P_hat`.
- `kappa = 1` makes `d^2 - v` unbiased for `delta^2`; larger `kappa` shrinks harder toward 0.
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
(cell 7: +.013 vs -.108), worse in the regular cells 2 and 4 (+.074 and +.028 vs +.004 and
-.001). Stronger shrinkage lowers the zero-contrast coefficient by at most about 2x and the
cell 5 and 6 bias by at most about 1.7x, while the regular-cell bias grows about tenfold.

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
