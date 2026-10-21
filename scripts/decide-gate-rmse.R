#!/usr/bin/env Rscript
# Decision rule for an RMSE gate with Monte Carlo error (Codex review G2, finding 5). PROPOSED, not registered.
# Usage: Rscript scripts/decide-gate-rmse.R <rmse> <mcse> <tolerance>
#   PASS     rmse + 2 mcse <= tolerance        clearly inside
#   FAIL     rmse - 2 mcse >  tolerance        clearly outside
#   CONFIRM  otherwise                         within 2 MCSE of the tolerance: run the independent confirmation block;
#                                              the verdict is decided by the confirmation run alone (same rule applied once).
# MCSE of an RMSE from R replicates with squared errors e2: sd(e2) / (2 * rmse * sqrt(R)) (delta method).
decide <- function(rmse, mcse, tol) {
  if (rmse + 2 * mcse <= tol) "PASS" else if (rmse - 2 * mcse > tol) "FAIL" else "CONFIRM"
}
mcse_rmse <- function(e, R = length(e)) { e2 <- e^2; stats::sd(e2) / (2 * sqrt(mean(e2)) * sqrt(R)) }
if (sys.nframe() == 0L) {
  # controls: planted clear pass, clear fail and borderline; and the MCSE formula against a Gaussian error sample
  stopifnot(decide(.10, .003, .15) == "PASS", decide(.20, .003, .15) == "FAIL", decide(.148, .003, .15) == "CONFIRM")
  set.seed(1); e <- stats::rnorm(2000, 0, .12)
  chk <- stats::sd(replicate(500, { x <- stats::rnorm(2000, 0, .12); sqrt(mean(x^2)) }))
  stopifnot(abs(mcse_rmse(e) / chk - 1) < .15)
  args <- commandArgs(TRUE)
  if (length(args) == 3L) cat(decide(as.numeric(args[1]), as.numeric(args[2]), as.numeric(args[3])), "\n")
  else cat("controls pass (pass / fail / confirm planted; MCSE formula within 15% of replicated SD)\n")
}
