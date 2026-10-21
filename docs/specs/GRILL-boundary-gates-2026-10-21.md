# GRILL: boundary point-risk plan and the per-cell gates to pre-register (2026-10-21)

Target: `tasks/plan-boundary-point-risk.md`, `docs/specs/SPEC-boundary-candidate-T3-2026-10-21.md`.
Interrogation of the author. Decisions below are the author's; none changes a frozen gate
until the author pre-registers it. Ledger appended after each branch.

## Decisions

### Q1. What anchors the new gates so they cannot be read as retrofitted?
**Decision: derived minimax floor.** Per-cell risk gates are set as a fixed multiple of a
two-point lower bound computed from the DGP laws alone (no candidate, no archived estimates).
Rejected: non-inferiority to k4 (k4 already fails the old gates; the margin is arbitrary),
dropping cells 5 and 7 (reads as relaxation), interval-only (leaves the point estimator open).
Open consequence: the floor and the multiple must be derived and stated before any fresh data.

### Q2. What makes per-cell risk gates meaningful given the Hodges problem?
**Decision: weak-signal ladder.** Add cells with the direct effect at 0 and about .25, .5 and 1
standard errors (fresh DGPs) and gate the maximum risk over the ladder against the floor.
Rejected: paired existing cells (spacing .23 SE is the spacing that defeats the old gate, so it
may be too tight and guards only two neighborhoods), per-cell-only with a documented loophole
(would pass a degenerate snap-to-boundary estimator).
Open consequence: ladder positions (in SE units, not raw coefficients, at fixed n) and the zero
indirect side (cell 6 analog) must be fixed before data; cost is a few extra cells in T6.

### Q3. Who picks the multiple, and from what?
**Decision: scientific tolerance, then feasibility check.** The author states the largest
acceptable error in P (a maximum RMSE over the ladder) on substantive grounds, without
reference to any candidate's archived risk. It is then checked against the derived floor and
registered. Rejected: an a priori multiple of the floor (still a judgment number, and the floor
depends on ladder design), calibrating on a disjoint pilot (fits the gate to the candidates).
Open consequence: the author must commit to the number; if it is below the floor the result is
a finding, not a reason to retune.

### Q4. What are PRs #54/#55 judged against once new gates exist?
**Decision: original gates stand for existing PR content.** The new gates apply only to a new
point estimator and its own change. The existing k4/selection work in #54/#55 stays judged by
the unchanged gates and stays held. Rejected: re-judging the PRs under the new gates (the new
gates would visibly decide a pre-existing result, the retrofit Q1 guards against); separating
the interval from the point estimator (changes the meaning of "complete diagnostics"; needs
the author's separate, explicit decision and a code-level split).
Open consequence: the PRs may stay blocked unless a new estimator, or a separate decision on
the interval, is made.

### Q5. How does the tolerance scale with n?
**Decision: rate-scaled tolerance per n.** Register the tolerance at n = 500 and scale it as
n^(-1/4) for n = 1000 (the D1 rate; local slopes .445 to .479 toward .5 in T2b), recomputing the
floor at each n. Rejected: n = 500 only (cells at n = 1000 get no gate), fixed absolute tolerance
across n (passes more easily at large n for reasons unrelated to the estimator).
Open consequence: commits to the D1 rate as an assumption; D3 (no threshold improves the rate)
is still an argument, not a theorem.

### Q6. What does the new gate set claim about the SE?
**Decision: report the SE only away from the boundary.** Keep SE/SD in [.90, 1.10] on the
regular cells (archived ratios about 1.03 to 1.07), and require the SE to be withheld (NA) when
a pre-registered boundary flag is on, with a gate on the flag's miss rate. Rejected: dropping the
SE from the point estimator (loses a calibrated SE where it works), keeping SE/SD in every cell
(feasibility of a calibrated SE near the boundary is unproven; likely a second infeasible gate).
Open consequence: the flag is itself a data-dependent threshold and needs its own definition
and validation; it can fail silently near the boundary, which is the Hodges problem again.

## Open Questions (not yet grilled)
- Ladder positions in SE units, and the zero-indirect side (cell 6 analog).
- The floor derivation per cell and how it is computed from the DGP laws (the TV method of
  `scripts/check-gate-feasibility-tv.R` generalizes to each ladder pair).
- Whether bias/SD, SE/SD and coverage gates are kept alongside the risk gate or replaced
  (SE calibration was never shown feasible near the boundary; coverage gate passes for Fieller).
- The boundary flag: its definition (threshold in SE units), tolerable miss rate on the ladder, and how its own Hodges-type failure is guarded.
- Review B timing: Codex review of the new gate text before any freeze.
- Hopper vs local for T6.
