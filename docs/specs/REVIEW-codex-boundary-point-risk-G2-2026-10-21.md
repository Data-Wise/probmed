# Codex adversarial review G2: complete gate set v1 (2026-10-21)

Target: `feature/boundary-point-risk` against `origin/dev`, focus `SPEC-boundary-gate-set-2026-10-21.md`. Review only; Codex applied nothing.

## Output (verbatim)

```text
# Codex Adversarial Review

Target: branch diff against origin/dev
Verdict: needs-attention

Do not register G4 yet. The proposed gates leave interval and SE failures untested, and the floor does not establish that the RMSE tolerances are feasible.

Findings:
- [high] Coverage gate names the wrong interval level (docs/specs/SPEC-boundary-gate-set-2026-10-21.md:48)
  Gate 6 calls for 90% interval coverage in [.92, .98]. The referenced validation script constructs 95% intervals with qnorm(0.975) and applies that range. If the new study constructs 90% intervals as written, a correctly calibrated method is expected to fail.
  Recommendation: Specify 95% intervals, or justify and set a coverage range appropriate to 90% intervals before registration.
- [high] Flagged ladder fits escape SE calibration (docs/specs/SPEC-boundary-gate-set-2026-10-21.md:61-65)
  The spec expects the flag to withhold SEs in roughly 88%–95% of ladder fits, yet gates SE/SD only at two regular cells and flag misses only at exact boundaries. An unreliable SE can therefore be exposed in unflagged interior fits without failing a gate. The Gaussian tables and n = 1000 SE scaling do not measure that conditional calibration.
  Recommendation: Predefine and gate SE calibration and flag behavior conditional on unflagged ladder fits at both sample sizes, using fresh DGP results.
- [high] The floor does not establish the ROOM verdicts (docs/specs/SPEC-boundary-gate-set-2026-10-21.md:44-47)
  Le Cam’s two-point floor is a necessary lower bound, not an attainable risk. The checker defines ROOM at an arbitrary 1.5 times that bound; the n = 1000 tolerances are only 1.52 times their reported floors. Thus the four ROOM labels do not establish that any estimator can meet the gates, while the pairwise calculation itself acknowledges that a stronger multipoint bound is uncomputed.
  Recommendation: Describe ROOM solely as “not ruled out,” and obtain an independent upper bound or reference-estimator feasibility check before presenting these tolerances as feasible.
- [medium] Nine grid points do not control risk between them (docs/specs/SPEC-boundary-gate-set-2026-10-21.md:52-53)
  The claim that the RMSE gate controls risk between nodes is unsupported. It tests the listed midpoints, but an estimator could have high risk elsewhere between those points; neither the two-point floor nor smoothness without a quantitative bound excludes that failure.
  Recommendation: Limit the claim to the nine tested points, or register a quantitative between-point risk bound and a check of its assumptions.
- [medium] Monte Carlo error can determine a close pass (docs/specs/SPEC-boundary-gate-set-2026-10-21.md:52)
  The spec acknowledges roughly 1.6%–2% relative RMSE MCSE at R = 2000, yet uses the observed RMSE as a hard threshold with no decision rule for close results. An estimator whose true risk is near a tolerance could pass or fail because of simulation noise.
  Recommendation: Before registration, define an independent confirmation run or a confidence-bound decision rule for results near a tolerance.

Next steps:
- Owner: author. Resolve these gate definitions and evidence gaps before Checkpoint G; retain the proposed status until the revised gate set is reviewed.
```

## Triage (proposed 2026-10-21; edits to the gate set await the author)

| # | Severity | Finding | Verdict | Proposed action |
|---|----------|---------|---------|-----------------|
| 1 | high | Gate 6 says 90% intervals | Accept, verified | `validate-wasserstein-boundary-selection.R` builds intervals with `qnorm(.975)` (lines 47-48, 57-58), so the original coverage gate is for 95% intervals. The "90%" in the SPEC is my error. Change to 95% intervals, range [.92, .98]. |
| 2 | high | SE/coverage not gated on unflagged ladder fits | Accept | Add a reported-and-gated conditional check: SE/SD in [.90, 1.10] and coverage in [.92, .98] among UNFLAGGED fits at ladder points with at least a stated number of unflagged fits (e.g. 200); the author sets the minimum and whether it is gated or only reported. |
| 3 | high | Floor does not establish ROOM | Accept | Reword ROOM as "not ruled out" everywhere (SPEC, checker output, todo). Add a feasibility reference: an independent upper-bound check, for example the achievable RMSE of the plug-in and of one regularized reference estimator on the DGP (no candidate tuning), before presenting tolerances as feasible. The multipoint bound stays an open limit. |
| 4 | medium | Grid does not control risk between nodes | Accept | Limit the claim to the nine tested points per path; optionally add a quantitative smoothness check (RMSE at extra random in-between points, reported). |
| 5 | medium | MC error decides close passes | Accept | Add a decision rule before registration: a pass needs observed RMSE <= tolerance; a result within 2 MCSE of a tolerance triggers an independent confirmation run on a disjoint seed block, decided by the confirmation alone. |
