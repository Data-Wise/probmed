# Codex adversarial review A: boundary point-risk (2026-10-21)

Plan: `tasks/plan-boundary-point-risk.md`. Target: working tree of `feature/boundary-point-risk` (7 untracked files). Output below is verbatim; triage follows it.
Codex could read the archived `replicates.csv` at the PR #55 checkout path given in the focus text. It saw no fresh data (none exists).

## Codex output (verbatim)

Target: working tree diff
Verdict: needs-attention

Brief: Review A blocks Checkpoint A. The threshold claim is false, and T2 does not independently test the claimed rate or risk trade-off.

Findings:
- [high] The weak-direct snap claim is false for λ = 1 (docs/specs/DERIVE-boundary-point-risk-2026-10-21.md:44-46)
  With the stated μ₂ = 0.23, D2's normal formula gives Pr(|Z + μ₂| ≤ 1) ≈ 0.670, not above 0.94. The archived 0.95 and 0.997 rates concern the larger tested thresholds. This overstates how broadly the cell 7 result applies.
  Recommendation: Restrict the claim to the tested thresholds and calculate snap probabilities from the stated formula for any other λ.
- [high] T2 simulates the estimator but does not test the derivation (scripts/check-boundary-risk-derivation.R:45-79)
  The prediction applies a copied point rule to Gaussian draws parameterized by archived truths and mean archived covariance. It never evaluates D1's coefficient or its n^(-1/4) rate across sample sizes, D2's independence-product snap formula, or a derived point-risk expression. Its bias and SD agreement therefore tests the Gaussian approximation and rule reproduction, not the central derivation; the SD gate itself fails 3/20.
  Recommendation: Add separately computed formula predictions, a sample-size or covariance-scaling check of the rate, and an explicit squared-risk decomposition with conditional snap terms. Keep the failed Gaussian check labeled as such.
- [medium] The cell 4 SD mechanism is claimed more precisely than the ablations establish (docs/specs/DERIVE-boundary-point-risk-2026-10-21.md:137-144)
  G1 and G4 reuse the archived d̂ distribution that produced the observed SD, while G3 changes covariance but only barely enters the gate. The script does not isolate tail shape, scale, and mean offset. Calling the cause identified risks carrying an unverified mechanism into T3's uncertainty design.
  Recommendation: Describe tail-distribution mismatch as a supported hypothesis; isolate its components with predictions specified independently of the observed cell 4 estimates before marking the mechanism identified.

Next steps:
- Next step — plan owner: correct D3 and revise Checkpoint A to reflect what T2 actually tests before proceeding to T3.

## Triage (2026-10-21)

| # | Severity | Finding | Decision | Action / reason |
|---|---|---|---|---|
| 1 | high | Weak-direct snap claim false for lambda = 1 | **Accept, fixed** | Verified: `Phi(.77) - Phi(-1.23) = .670` (also .944 at 1.96, .994 at 2.807). D3 in `DERIVE-...md` restricted to lambda of about 1.96 and above, with the formula stated. |
| 2 | high | T2 simulates the estimator, does not test the derivation (no rate check, no formula-vs-simulation check of D2, no risk decomposition) | **Accept, open** | True: the scripts compare a copied rule on Gaussian/archived draws. New task T2b added to `tasks/todo-boundary-point-risk.md`. DERIVE doc now states what T2 does and does not test. Unresolved until T2b runs, so it blocks Checkpoint A. |
| 3 | medium | Cell 4 mechanism overstated as "identified" | **Accept, fixed** | Reworded to "supported hypothesis"; scale/skew/mean offset not isolated; G1/G4 reuse the observed sample. Isolation deferred (not needed to proceed unless T3 relies on it; T3 must not). |

Unresolved high-severity findings: **1** (finding 2). Checkpoint A stays blocked.
No fix was applied by Codex; all edits above were made by us.

## Follow-up: finding 2 addressed by T2b (2026-10-21)

T2b (`scripts/check-boundary-risk-rate.R`) now tests the derivation directly: rate slope
.459 and leading coefficient .964 (pass), snap formula 10/10 vs exact and simulation, a
discriminating control, and a branch-wise RMSE decomposition that passes 8/10 and fails in
the two known rare-snap cells (cell 4 select95, cell 5 select995). Finding 2 is therefore
addressed, not "fully passed": C1 fails and D3 remains an argument. Whether this clears
Checkpoint A is the author's decision.
