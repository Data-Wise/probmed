# G4: complete gate set v1 for a new binary Wasserstein point estimator (2026-10-21)

**Status: PROPOSED for Codex review G2. NOT REGISTERED, NOT FROZEN.** Registration is Checkpoint G (a commit before any
fresh data). Entries marked "author" were entered by the author on 2026-10-21 by adopting the assistant's recommendation
and are revisable until then; entries marked "assistant" are drafting choices the author must confirm at Checkpoint G.
Inputs: `GRILL-boundary-gates-2026-10-21.md`, `REVIEW-codex-boundary-point-risk-G-2026-10-21.md` and `-H-`,
`DERIVE-boundary-floor-2026-10-21.md`, `SPEC-boundary-candidate-T3-2026-10-21.md`.

**Scope.** Applies only to a NEW point estimator and its own change. PRs #54/#55 stay judged by their original frozen gates
(Q4) and stay held. No frozen diagnostic is changed. Estimand: binary-outcome Wasserstein proportion mediated
P = sqrt|d1| / (sqrt|d1| + sqrt|d2|), DGP of `scripts/validate-wasserstein-joint-contrast.R` (PR #55 branch).

## 1. Cells

| Ids | Cells | n |
|-----|-------|---|
| 1 | binary regular (original cell 2: alpha .7, direct .6, interaction .2) | 500 |
| 2 | binary regular (original cell 4) | 1000 |
| 3-11 | direct ladder, 9 points (alpha .7, interaction 0; direct = boundary 0 and the points below) | 500 |
| 12-20 | indirect ladder, 9 points (direct .6, interaction .2; alpha = boundary 0 and the points below) | 500 |
| 21-29 | direct ladder | 1000 |
| 30-38 | indirect ladder | 1000 |

Ladder points (`docs/specs/data/gate-grid-g1c-2026-10-21.txt`: boundary, four TV-matched nodes, midpoints of consecutive points):

| n | Path | Parameter values |
|---|------|------------------|
| 500 | direct | 0, .00875, .0175, .02635, .0352, .05325, .0713, .1095, .1477 |
| 500 | indirect | 0, .0040, .0080, .0120, .0160, .0242, .0324, .05005, .0677 |
| 1000 | direct | 0, .00625, .0125, .01875, .025, .03775, .0505, .0777, .1049 |
| 1000 | indirect | 0, .00285, .0057, .00855, .0114, .01725, .0231, .03525, .0474 |

Positions are set by data-law TV to the boundary law (.05, .10, .20, .40 at the nodes), not by SE units of d2 (audit finding).
R = 2000 replicates per cell, as in the original studies: assistant; confirm.
Seeds (assistant, proposed, not frozen): `seed = 1,000,000,000 + cell id * 10000 + replicate`, cell 1..38, replicate 1..2000;
audited disjoint from 13 earlier ranges by `scripts/audit-seed-ranges.R` (block [1000010001, 1000405000], allowing cell ids up to 40
and replicates up to 5000); re-audit at T5.

## 2. Gates

| # | Gate | Value | Entered by | Feasibility check |
|---|------|-------|------------|-------------------|
| 1 | zero failed fits, every cell | 0 | unchanged original gate | n/a |
| 2 | max RMSE over the direct ladder (9 points), n = 500 | <= .15 | author | G3 ROOM (floor .0926, 1.62x) |
| 3 | max RMSE over the indirect ladder, n = 500 | <= .12 | author | G3 ROOM (floor .0675, 1.78x) |
| 4 | max RMSE over the direct ladder, n = 1000 | <= .125 | author | G3 ROOM (floor .0820, 1.52x) |
| 5 | max RMSE over the indirect ladder, n = 1000 | <= .09 | author | G3 ROOM (floor .0591, 1.52x) |
| 6 | 90% interval coverage in [.92, .98] at the regular cells (ids 1-2) and the two boundary points per n (direct 0, alpha 0) | unchanged original range | assistant (scope) | n/a (Fieller passes coverage in the archived binary cells) |
| 7 | mean SE / empirical SD in [.90, 1.10], regular cells (ids 1-2) only | unchanged original range | author (Q6, scoped to regular cells) | n/a |
| 8 | boundary flag miss rate at the two boundary points per n <= .06 | see section 3 | author | tables in `tabulate-flag-misses-2026-10-21.txt` |

No rounding at any gate. The RMSE gates use the observed Monte Carlo RMSE with its MCSE reported alongside; no allowance is added.
RMSE is over all 9 grid points per path at that n, so the gate also controls risk between the nodes (G1c).

**Dropped as a gate for the new estimator:** |bias| / empirical SD < .10. It is infeasible in cells 5 and 7 (TV .1404 vs needed .179;
T3 SPEC section 2) and can be gamed by variance inflation outside [0, 1]. It is REPORTED at every cell, not gated.
**Reported, not gated:** SE/SD and coverage at interior ladder points, the flag OFF rate at every grid point, selection or snap rates.

## 3. Boundary flag (gate 8)

ON when min over j of |d_hat_j| / se_j < 1.96. At the two boundary points (direct 0, alpha 0) the flag stays OFF (a miss) in
at most 6% of replicates (5% nominal + 2 MCSE at R = 2000; the archive shows .056 and .055). The SE is withheld (NA) when the
flag is ON; the SE/SD gate (gate 7) applies only to the regular cells. The flag does NOT protect any ladder point: at c = 1.96 it is
ON for roughly 88% to 95% of ladder fits, so the SE is withheld there by design, and ladder risk is gated by RMSE alone.
Carried-forward caveats: Gaussian-approximation tables, n = 1000 SE scaled by 1/sqrt(2), indirect-path SE not constant (.0172 to .0235).

## 4. Procedure

1. Review G2 (Codex) on this document before registration; triage every finding.
2. Checkpoint G: the author registers this table in a commit before any fresh data; feasibility verdicts are re-run if any
   tolerance changes.
3. T4 tests, Review B, T5 freeze (re-audit the seed block, hashes), then T6 on the fresh seeds (38 cells x 2000 replicates; the
   earlier study ran 44,000 rows in 318 s locally, so local vs hopper is a T6 decision).
4. A failed feasibility check is a finding, not a reason to retune a tolerance. If a candidate cannot meet the gates, or the gates
   are NEAR-FLOOR or INFEASIBLE, the endpoint is the negative result of the T3 SPEC section 5, option 3.

## 5. Ledger effect (applies at Checkpoint G)

Q1 (floor anchor): used as a feasibility check on author-set tolerances, not as a multiple. Q2 (ladder): implemented with TV-matched positions.
Q3 (tolerance first, then feasibility): followed (G3). Q4 (PRs #54/#55 under original gates): unchanged. Q5: reversed, per-n tolerances.
Q6 (SE only away from the boundary): scoped to regular cells with the flag of section 3. Promotion from PROVISIONAL happens at Checkpoint G.

## 6. Open items for review

- Coverage gate 6 scope (assistant): which cells carry the original coverage gate for the new estimator.
- R = 2000 per cell gives RMSE relative MCSE of roughly 1.6% to 2%; an RMSE within about 3% of a tolerance is a coin flip.
- The Gaussian-approximation basis of gate 8 and the n = 1000 SE scaling are unverified assumptions.
