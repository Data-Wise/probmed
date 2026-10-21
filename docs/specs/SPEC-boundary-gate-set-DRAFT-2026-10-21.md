# G4 DRAFT: one complete gate set (skeleton, 2026-10-21)

**Status: DRAFT. Gates 1-3 entered by the author 2026-10-21 (G1a, G1b) and pre-checked by G3: all ROOM (`docs/specs/data/check-gate-tolerance-2026-10-21.txt`). Not frozen; registration is Checkpoint G, after G1c, G1d, G4 and Review G2.** Every value marked AUTHOR is yours to write (G1). Nothing here is
frozen, nothing changes a frozen diagnostic, and PRs #54/#55 stay judged by the original gates (Q4).
Inputs: `GRILL-boundary-gates-2026-10-21.md`, `REVIEW-codex-boundary-point-risk-G-2026-10-21.md`,
`DERIVE-boundary-floor-2026-10-21.md`. Applies only to a NEW point estimator and its own change.

## 1. Estimand and cells

Binary-outcome Wasserstein proportion mediated P = sqrt|d1| / (sqrt|d1| + sqrt|d2|), same DGP as
`scripts/validate-wasserstein-joint-contrast.R` (PR #55 branch). Seeds disjoint from every earlier
study: PROPOSED (not frozen) seed = 1,000,000,000 + cell * 10000 + replicate, cell 1..20, replicate 1..5000, i.e. [1000010001, 1000205000]; `scripts/audit-seed-ranges.R` finds no overlap with 13 earlier ranges and no seed literal inside it (output `docs/specs/data/audit-seed-ranges-2026-10-21.txt`). Re-audit at T5 and freeze there.

## 2. Ladder (positions by data-law TV to the boundary law, from G2)

| Ladder | n = 500 positions at TV .05 / .10 / .20 / .40 | n = 1000 |
|--------|-----------------------------------------------|----------|
| direct (alpha .7, interaction 0) | .0175 / .0352 / .0713 / .1477 | .0125 / .0250 / .0505 / .1049 |
| indirect (direct .6, interaction .2) | .0080 / .0160 / .0324 / .0677 | .0057 / .0114 / .0231 / .0474 |

Plus the boundary points (direct 0, alpha 0) and the existing regular cells. Whether to keep all four
TV levels: AUTHOR. Ladder-to-ladder pairs are not yet in the floor (open).

## 3. Gates (one row per gate; Review G finding in brackets)

| # | Gate | Value | Basis |
|---|------|-------|-------|
| 1 | max RMSE over the direct ladder, n = 500 | **.15** (author, 2026-10-21) | substantive tolerance; G3 verdict must not be INFEASIBLE [Q1-Q3] |
| 2 | max RMSE over the indirect ladder, n = 500 | **.12** (author, 2026-10-21) | same |
| 3 | same, n = 1000 (own tolerances, not scaled) | **.125 direct, .09 indirect** (author, 2026-10-21) | Q5 reversed; floors .082 / .059 |
| 4 | max RMSE over the finer grid (boundary, 4 TV-matched nodes, midpoints of consecutive points; 9 per path, `docs/specs/data/gate-grid-g1c-2026-10-21.txt`) at the SAME tolerances as gates 1-3 | **chosen** (author, 2026-10-21, G1c) | control for estimators that pass only at the nodes [Review G]; the pairwise floors (DERIVE section 7) are not above the boundary-pair floor, so G3 verdicts stand |
| 5 | Fieller interval coverage | unchanged original gate | already passes in the binary cells |
| 6 | SE/SD on regular cells | [.90, 1.10] (Q6, provisional) | archived ratios about 1.03 to 1.07 |
| 7 | SE withheld (NA) when the boundary flag is on; flag miss rate | AUTHOR | Q6 [Review G] |

Dropped as a gate for the new estimator: bias/SD < .10 (infeasible in cells 5 and 7, gameable outside
[0,1]; T3 SPEC section 2). It may still be reported.

## 4. Boundary flag (needed by gate 7)

Definition (threshold in SE units), tolerable miss rate on the ladder, and how its own Hodges-type
failure is guarded: AUTHOR. Not specified here because any threshold I pick would be tuned.

## 5. Procedure constraints

- Freeze (T5) only after G3 returns no INFEASIBLE and Codex reviews this table; the freeze commit records
  the triage. Fresh seeds and archived rows per the plan; no gate edited after data exist.
- A failed feasibility check is a finding to report, not a reason to retune the tolerance.
- If every admissible tolerance is INFEASIBLE or NEAR-FLOOR, the endpoint is the negative result in
  the T3 SPEC section 5, option 3.
