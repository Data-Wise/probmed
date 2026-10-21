#!/usr/bin/env Rscript
# G3: feasibility pre-check of a proposed max-RMSE tolerance against the two-point floor.
# Usage: Rscript scripts/check-gate-tolerance.R <n> <tol_direct> [<tol_indirect>]
# Reads docs/specs/data/floor-ladder-n<n>.csv (boundary-to-node) AND floor-pairs-n<n>.csv (all node pairs), both
# written by scripts/check-boundary-floor-ladder.R and check-boundary-floor-pairs.R (Codex review 2026-10-21,
# finding 2: use the maximum supported floor across both). floor_hi is the same bound with TV lowered by 2 MC SE.
# Verdicts per path (tol is a bound on max RMSE):
#   INFEASIBLE   tol <= floor               no estimator can meet it (Le Cam), a finding, not a reason to retune
#   BORDERLINE   floor < tol <= floor_hi    within 2 MC SE of the floor; treat as INFEASIBLE unless a larger R says otherwise
#   NEAR-FLOOR   floor_hi < tol < 1.5 floor floor is a lower bound, not attainable; treat as likely infeasible
#   NOT-RULED-OUT tol >= 1.5 floor          the two-point floor does not exclude it; this is NOT evidence that any estimator can meet it (Codex review G2, finding 3)
# The 1.5 cut is a stated convention for reporting, not a derived constant.
# Positive controls run first: half the floor gives INFEASIBLE, twice the floor gives NOT-RULED-OUT, and a tolerance that
# sits below the pairwise maximum but above the ladder maximum gives INFEASIBLE (the case the review flagged).
args <- commandArgs(TRUE)
stopifnot(length(args) >= 2L)
n <- as.integer(args[1]); tol <- as.numeric(args[-1])
if (length(tol) == 1L) tol <- c(tol, tol)
rd <- function(f) {
  d <- utils::read.csv(sprintf("docs/specs/data/%s-n%d.csv", f, n))
  d$floor_hi <- d$gap * sqrt((1 - pmax(0, d$tv - 2 * d$se)) / 8)
  d[, c("path", "floor_rmse", "floor_hi")]
}
f <- rbind(rd("floor-ladder"), rd("floor-pairs"))
fl <- tapply(f$floor_rmse, f$path, max)
fh <- tapply(f$floor_hi, f$path, max)
verdict <- function(t, fl, fh) if (t <= fl) "INFEASIBLE" else if (t <= fh) "BORDERLINE" else if (t < 1.5 * fl) "NEAR-FLOOR" else "NOT-RULED-OUT"
lad <- tapply(rd("floor-ladder")$floor_rmse, rd("floor-ladder")$path, max)
pr <- tapply(rd("floor-pairs")$floor_rmse, rd("floor-pairs")$path, max)
mid <- (lad[["direct"]] + pr[["direct"]]) / 2
stopifnot(verdict(.5 * fl[["direct"]], fl[["direct"]], fh[["direct"]]) == "INFEASIBLE",
  verdict(2 * fl[["direct"]], fl[["direct"]], fh[["direct"]]) == "NOT-RULED-OUT",
  pr[["direct"]] < lad[["direct"]] || verdict(mid, fl[["direct"]], fh[["direct"]]) == "INFEASIBLE")
for (i in 1:2) {
  path <- c("direct", "indirect")[i]
  cat(sprintf("n = %d, %-8s floor %.4f (hi %.4f), tolerance %.4f, ratio %.2f: %s\n",
    n, path, fl[[path]], fh[[path]], tol[i], tol[i] / fl[[path]], verdict(tol[i], fl[[path]], fh[[path]])))
}
