#!/usr/bin/env Rscript
# G3: feasibility pre-check of a proposed max-RMSE tolerance against the two-point floor.
# Usage: Rscript scripts/check-gate-tolerance.R <n> <tol_direct> [<tol_indirect>]
# Reads docs/specs/data/floor-ladder-n<n>.csv written by scripts/check-boundary-floor-ladder.R.
# Verdicts per path (floor = largest two-point floor on that ladder; tol is a bound on max RMSE):
#   INFEASIBLE   tol <= floor            no estimator can meet it (Le Cam), a finding, not a reason to retune
#   NEAR-FLOOR   floor < tol < 1.5 floor floor is a lower bound, not attainable; treat as likely infeasible
#   ROOM         tol >= 1.5 floor        not excluded by the floor (says nothing about any candidate)
# The 1.5 cut is a stated convention for reporting, not a derived constant.
# Positive controls run first: a tolerance of half the floor must give INFEASIBLE, twice the floor ROOM.
args <- commandArgs(TRUE)
stopifnot(length(args) >= 2L)
n <- as.integer(args[1]); tol <- as.numeric(args[-1])
if (length(tol) == 1L) tol <- c(tol, tol)
f <- utils::read.csv(sprintf("docs/specs/data/floor-ladder-n%d.csv", n))
fl <- tapply(f$floor_rmse, f$path, max)
verdict <- function(t, fl) if (t <= fl) "INFEASIBLE" else if (t < 1.5 * fl) "NEAR-FLOOR" else "ROOM"
stopifnot(verdict(.5 * fl[["direct"]], fl[["direct"]]) == "INFEASIBLE",
  verdict(2 * fl[["direct"]], fl[["direct"]]) == "ROOM")
for (i in 1:2) {
  path <- c("direct", "indirect")[i]
  cat(sprintf("n = %d, %-8s floor %.4f, tolerance %.4f, ratio %.2f: %s\n",
    n, path, fl[[path]], tol[i], tol[i] / fl[[path]], verdict(tol[i], fl[[path]])))
}
