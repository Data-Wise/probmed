# Local CI transition (2026-10-20)

## Decision and live constraints

Use local package checks for PR validation and post the result as the required
`local-ci/probmed-full` commit status. On 2026-10-20, `main` still required
`ubuntu-latest (release)` with strict up-to-date checking. `dev` had no branch
protection, and there were no open PRs. Do not remove that old requirement until
the replacement script is integrated, exercised on a real PR, and its status is
visible for the PR head. Neither a local pass nor a status authorizes a merge.

## Gates

- `scripts/check-local.sh` remains the fast, non-authoritative preflight.
- `scripts/check-full-local.sh` requires a clean committed checkout and runs
  the `DESCRIPTION`-pinned roxygen2 version, `devtools::document()` without
  generated drift, tests with `NOT_CRAN=true`, `rcmdcheck` with
  `--no-manual`, `covr` (overall coverage >90%, no uncovered expressions in
  `R/compute-core.R`), and pkgdown in a disposable copy that excludes agent
  instruction files. Logs remain under the printed temporary directory.
  Use an isolated R library via `R_LIBS_USER` for the pinned roxygen2 7.3.3;
  the current system version 8.1.0 rewrites unrelated Rd files. The scoped
  Rd synchronization in this branch was generated with 7.3.3.
- For an open PR, run `scripts/check-full-local.sh --pr NUMBER` from its exact
  head. It checks the remote, PR base/head/state, rejects fork PRs, posts
  `pending`, and posts
  `success` only if every gate passes and the head has not moved. Failures post
  `failure`. The GitHub CLI must have permission to write commit statuses.
- This local Mac run does **not** reproduce the old Windows/Linux/R-version
  matrix or Codecov upload. Use the manual R-CMD-check, coverage, and R-hub
  workflows when cross-platform or release evidence is needed. Do not label
  the local status as a matrix pass.

## Safe switch order

1. Integrate this feature into `dev` through a reviewed PR. Keep the live
   `main` protection unchanged. Coverage becomes manual, while R-CMD-check
   keeps its `main` PR trigger until the status switch; pkgdown keeps
   `main`/release publication.
2. On a future `dev` to `main` PR, run the full local gate on that PR head and
   verify the `local-ci/probmed-full` status exists on the same SHA. The old
   Ubuntu check stays required and runnable during this transition.
3. Only then replace `main`'s required `ubuntu-latest (release)` context with
   `local-ci/probmed-full`, preserving strict mode. Verify the live protection
   response and PR mergeability. If any check is missing, restore the original
   context before proceeding.
4. In a follow-up reviewed change, remove R-CMD-check's `main` PR trigger,
   leaving manual dispatch. Review and merge only with separate author
   authorization. Never bypass a required check. Keep `pkgdown` deployment
   and manual R-hub available.

This is a reduction in independent platform coverage and automated execution.
The status is posted with a user's GitHub credentials, not by a trusted hosted
runner; branch protection should restrict who may write the status where
possible. Keep that tradeoff explicit in release and submission review.
