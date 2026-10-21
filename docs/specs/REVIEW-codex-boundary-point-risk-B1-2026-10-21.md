# Codex adversarial review B1: gate-infeasibility argument (2026-10-21)

Plan: `tasks/plan-boundary-point-risk.md`. Scope: `docs/specs/SPEC-boundary-candidate-T3-2026-10-21.md` section 4, `scripts/check-gate-feasibility-tv.R`, grid claims in section 3. Requested before the author proposes feasible gates. Output verbatim; triage follows. Codex could read the PR #55 validation report and DGP script at the checkout path given.

## Codex output (verbatim)

Target: working tree diff
Verdict: needs-attention

Do not ship the infeasibility result as a claim about the frozen empirical gate. The population TV argument is sound for a common [0,1]-valued estimator, including randomized estimators, but the stated consequence exceeds it.

Findings:
- [high] Population impossibility is presented as an empirical gate impossibility (docs/specs/SPEC-boundary-candidate-T3-2026-10-21.md:75-84)
  The frozen validator uses bias and sample SD from 2,000 replicates. Section 4 instead bounds population means and SDs. A finite sample can pass both empirical gates even when the population inequalities cannot hold; sample SD is not bounded by sqrt(E(1-E)). The claimed impossibility 'by construction' and decision to skip validation do not follow without a bound on empirical pass probability.
  Recommendation: State the result as population gate infeasibility and quantify the chance of an empirical pass before using it to rule out T6.
- [medium] Smooth candidate comparison contradicts its own grid (docs/specs/SPEC-boundary-candidate-T3-2026-10-21.md:50-56)
  The claim that the family is 'no better than k4 at the stress cells' conflicts with the displayed biases: at kappa=3.84, tau=.5, cell 7 bias is +.013 versus k4's -.108; cells 5 and 6 also improve in absolute bias. This can wrongly dismiss a candidate despite the grid's stated purpose of guiding selection.
  Recommendation: Compare the family with k4 cell by cell and separate bias improvements from failure of the complete gate.
- [medium] Claimed TV controls are absent from the check (docs/specs/SPEC-boundary-candidate-T3-2026-10-21.md:71-73)
  Section 4 reports TV controls for shifts 0, .25, and .5, but the script evaluates only the .05 shift and its saved output contains no controls. The reported controls therefore cannot be reproduced from the cited check, weakening the numerical evidence for the central no-go decision.
  Recommendation: Implement and record the three controls in the cited script, or remove the unsupported claim.

Next steps:
- Owner: plan author. Correct the empirical-gate conclusion and grid comparison, then rerun the TV check with its stated controls.

## Triage (2026-10-21)

| # | Severity | Finding | Decision | Action |
|---|---|---|---|---|
| 1 | high | Population impossibility presented as empirical-gate impossibility | **Accept, fixed** | The sample SD of [0,1] values satisfies `s^2 <= R/(R-1) m(1-m)`, so the same algebra holds for sample means; an empirical pass needs noise >= .039 (3.4 SD). Bernstein bound on P(empirical pass) <= .0042 (normal approx .0003). "By construction" withdrawn; spec section 4 and the recommendation reworded. Verdict "unattainable" is now "population-infeasible; empirical pass probability < .005". |
| 2 | medium | "No better than k4 at the stress cells" contradicts the grid | **Accept, fixed** | Verified against the table (kappa 3.84, tau .5: cell 7 +.013 vs -.108). Replaced with a cell-by-cell comparison. |
| 3 | medium | TV controls claimed but not in the script | **Accept, fixed** | Controls (shift 0, .25, .5) added to the script and its saved output (0, .623, .921). Main TV re-estimated at .1404 (earlier .1414, RNG stream). |

Codex's own assessment: the population TV argument is sound for a common [0,1]-valued
estimator, including randomized estimators. No unresolved high-severity finding. Codex
applied no fixes.
