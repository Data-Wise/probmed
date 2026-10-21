# Codex adversarial review H: gate tooling (G2, G2b, G2c, G3) (2026-10-21)

Target: `feature/boundary-point-risk` against `origin/dev` (39 files). Review only; Codex applied nothing.
Scope note from Codex: the T2 derivation check could not run in its read-only sandbox (no R temp dir); that is a sandbox limit, not a finding about the script.

## Output (verbatim)

```text
# Codex Adversarial Review

Target: branch diff against origin/dev
Verdict: needs-attention

Do not ship the gate tooling yet: two checks can report a reassuring result while omitting evidence already required by the plan. The T2 derivation check could not run because the read-only sandbox could not create an R temporary directory.

Findings:
- [medium] Seed audit can pass without the PR #55 scripts (scripts/audit-seed-ranges.R:35-50)
  The PR #55 directory is optional in the discovery scan. If its default checkout is absent or moved, the script skips it and can still print AUDIT PASS, so the proposed fresh block has not been checked against that source.
  Recommendation: Require the PR #55 directory to exist and contain the expected scripts; fail the audit when that input is unavailable.
- [medium] Tolerance check ignores the newer pairwise floors (scripts/check-gate-tolerance.R:15-17)
  G3 reads only floor-ladder CSVs. The committed n=500 pairwise results have larger floors: 0.09257 versus 0.09226 direct and 0.06747 versus 0.06596 indirect. A tolerance between either pair of values would be reported as NEAR-FLOOR although the branch's own pairwise estimate marks it INFEASIBLE.
  Recommendation: Use the maximum supported floor across ladder and pairwise results for each path and sample size, and account for Monte Carlo uncertainty near the cutoff.

Next steps:
- Owner: branch author. Fix both gate checks before using them to freeze tolerances or seeds.
```

## Triage

| # | Finding | Decision | Action and evidence |
|---|---------|----------|---------------------|
| 1 | Seed audit passes without PR #55 scripts | Accept | `audit-seed-ranges.R` now exits 2 with AUDIT FAIL unless the directory and the four expected scripts exist. Control: `/nonexistent` gives rc 2; the real directory still gives AUDIT PASS. |
| 2 | G3 ignores pairwise floors; no MC allowance | Accept | `check-gate-tolerance.R` takes the max floor over `floor-ladder` and `floor-pairs` CSVs and adds a BORDERLINE verdict within 2 MC SE of the cutoff (`floor_hi`). Control: a tolerance between the ladder and pairwise maxima (.0925 direct) now reads INFEASIBLE. n = 500 floors are now .0926 / .0675. |

Side effect on advice: at n = 1000 the direct floor is .0820, so tolerance .12 is 1.46x (NEAR-FLOOR); .125 or more clears the 1.5x convention.
