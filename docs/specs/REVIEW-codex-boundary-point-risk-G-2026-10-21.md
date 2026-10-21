# Codex adversarial review G: gate-redesign grill (2026-10-21)

Target: `docs/specs/GRILL-boundary-gates-2026-10-21.md` (Q1-Q6), run with `--base ecc982d` (the grill commit `e070948` only). Output verbatim; triage follows. The grill decisions are the author's; this record does not reverse any of them.

## Codex output (verbatim)

Target: branch diff against ecc982d
Verdict: needs-attention

Do not freeze these gates. Reverse Q1–Q3 as the proposed acceptance rule, Q5, and Q6. Keep Q4’s hold on PRs #54/#55.

Findings:
- [high] The proposed floor does not establish a non-retrofitted, feasible gate (docs/specs/GRILL-boundary-gates-2026-10-21.md:9-32)
  The two-point result in T3 is an infeasibility bound for the old bias gate, not a minimax risk calculation for the proposed ladder. Q1 leaves both the risk floor and its multiplier uncomputed; Q3 then chooses a tolerance after the archived candidate biases and troublesome cells are known. Excluding archived estimates from the floor calculation does not remove that specification choice. A finite ladder also cannot, by itself, establish risk control between its points. The resulting pass rule could be tailored to known failures while being described as independently derived.
  Recommendation: Reverse Q1–Q3 as an acceptance rule. Derive the stated risk bound, justify the tolerance independently, and specify how risk between ladder points is bounded before freezing a gate.
- [high] The sample-size scaling applies a raw plug-in bias rate to a different risk target (docs/specs/GRILL-boundary-gates-2026-10-21.md:44-50)
  D1 supports an asymptotic n^(-1/4) bias at an exact zero contrast for the raw plug-in. Q5 applies that rate to maximum RMSE for a new estimator across a weak-signal ladder at n = 500 and 1000. Neither D1 nor the reported local slopes establishes that finite-sample RMSE scaling; D3 remains an argument. This can make the n = 1000 pass threshold too permissive or too strict without revealing which.
  Recommendation: Reverse Q5. Set separate, substantively justified tolerances at each n, or derive and independently validate the proposed RMSE scaling for the estimator and ladder being gated.
- [high] The boundary flag can make SE calibration selectively untested (docs/specs/GRILL-boundary-gates-2026-10-21.md:52-59)
  Q6 checks SE/SD only in regular cells and requires NA when a data-dependent flag fires. Its flag definition and miss-rate tolerance remain open, with no stated calibration check for unflagged weak-signal cases or limit on withholding. An estimator could suppress SEs in difficult cases while passing the regular-cell ratio, and a flag miss could expose an unreliable SE to users.
  Recommendation: Reverse Q6 as a gate decision until the flag, false-negative and withholding limits, and SE calibration conditional on an unflagged result are specified and tested across the ladder.
- [medium] The freeze procedure still forbids the proposed gate change (docs/specs/GRILL-boundary-gates-2026-10-21.md:61-68)
  The plan requires no gate relaxation, an unchanged-gate T6 run, and a verdict under every unchanged gate. This ledger proposes new gates but leaves Review B timing and whether existing bias/SD, SE/SD, and coverage gates remain open. A T5 freeze under the current plan would therefore have no unambiguous pass rule. Fresh seeds prevent reusing datasets; they do not resolve this already informed choice of diagnostics.
  Recommendation: Before T5, amend the plan and spec with one complete gate set, its scope, and Review B approval of that exact text. Keep PRs #54/#55 under their original gates as Q4 states.

Next steps:
- Owner: author. Decide whether to pursue a newly justified point-risk gate, then reconcile the plan and spec before any freeze or fresh-seed run.

## Triage (2026-10-21)

Codex applied no fixes. Status below is my assessment; each reversal is the author's call.

| # | Sev | Finding | Assessment | Proposed handling |
|---|---|---|---|---|
| 1 | high | Q1-Q3: floor uncomputed, tolerance set after the failures are known; a finite ladder does not bound risk between points | **Valid.** The floor and multiple are not yet derived, so "non-retrofitted" is overstated; the author's tolerance will be an informed choice. My Q2 claim that the ladder "closes the loophole" was too strong. | Mark Q1-Q3 provisional until (a) the floor is derived for the actual ladder, (b) the tolerance is justified independently of the archived risks and the informed-choice caveat is disclosed, (c) between-point risk is bounded (denser ladder or a smoothness argument). |
| 2 | high | Q5: n^(-1/4) bias rate applied to maximum RMSE of a new estimator | **Valid.** D1 is a bias rate for the raw plug-in at an exact zero; D3 is an argument. My recommendation was weak here. | Reverse Q5 as written: separate, substantively justified tolerances at n = 500 and 1000, unless a RMSE rate is derived and validated. |
| 3 | high | Q6: SE withheld near a data-dependent flag can hide failures | **Valid.** Flag, false-negative limit, withholding cap and SE calibration conditional on unflagged are all unspecified. | Mark Q6 provisional until those four are specified and tested across the ladder. |
| 4 | med | The plan still requires unchanged gates and no relaxation; no complete new gate set exists | **Valid.** T5 under the current plan would have no unambiguous pass rule. | Before T5, amend the plan and T3 spec with ONE complete gate set and scope, and run Review B on that exact text. Q4 (PRs stay under original gates) stands, as Codex agrees. |

Applied 2026-10-21 on the author's instruction: Q5 reversed; Q1-Q3 and Q6 marked provisional in the ledger.
