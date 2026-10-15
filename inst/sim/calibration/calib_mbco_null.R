# Composite-null calibration harness for method = "mbco" (Phase 8, T2).
#
# Usage (from the package root):
#   Rscript inst/sim/calibration/calib_mbco_null.R <task_id> <outdir> [reps] [width_reps]
# Lattice: (a, b) = (h_a, h_b) / sqrt(N), h in {0, 0.5, 1, 2, 4}^2, N in {100, 1000},
# covariate in {FALSE, TRUE} -> 100 points; task_id 1..100 runs one point.
#
# Coverage of an LR-inversion interval equals acceptance of the truth, so each
# replication records LR at the true P_med plus the two branch statistics
# (H0: a = 0, H0: b = 0). Any cutoff rule c(LR_a0, LR_b0) is then scored
# offline as mean(LR_true <= c). The first `width_reps` replications also run
# the full interval, as a self-check: interval coverage must equal acceptance
# coverage at the chi-square(1) cutoff.
#
# Set PROBMED_LOAD=source to use pkgload::load_all(".").

args <- commandArgs(trailingOnly = TRUE)
task <- as.integer(args[1]); outdir <- args[2]
reps <- if (length(args) >= 3) as.integer(args[3]) else 2000L
width_reps <- if (length(args) >= 4) as.integer(args[4]) else 50L
if (identical(Sys.getenv("PROBMED_LOAD"), "source")) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else suppressMessages(library(probmed))
ns <- asNamespace("probmed")

grid <- expand.grid(h_a = c(0, 0.5, 1, 2, 4), h_b = c(0, 0.5, 1, 2, 4),
                    N = c(100L, 1000L), cov = c(FALSE, TRUE))
pt <- grid[task, ]
a <- pt$h_a / sqrt(pt$N); b <- pt$h_b / sqrt(pt$N)
truth <- stats::pnorm(a * b / sqrt(2 * (b^2 + 1)))   # Vm = Vy = 1, x contrast 1

RNGkind("L'Ecuyer-CMRG"); set.seed(20261015L)
s <- .Random.seed
for (i in seq_len(task)) s <- parallel::nextRNGStream(s)

one <- function(r) {
  d <- data.frame(X = stats::rbinom(pt$N, 1, 0.5), C = stats::rnorm(pt$N))
  d$M <- a * d$X + (if (pt$cov) 0.5 * d$C else 0) + stats::rnorm(pt$N)
  d$Y <- 0.3 * d$X + b * d$M + (if (pt$cov) 0.3 * d$C else 0) + stats::rnorm(pt$N)
  fm <- if (pt$cov) M ~ X + C else M ~ X
  fy <- if (pt$cov) Y ~ X + M + C else Y ~ X + M
  ex <- medfit::extract_mediation(stats::lm(fm, d), model_y = stats::lm(fy, d),
                                  treatment = "X", mediator = "M", data = d)
  prep <- ns$.mbco_prep(ex, 0, 1)
  out <- data.frame(task, h_a = pt$h_a, h_b = pt$h_b, N = pt$N, cov = pt$cov,
                    rep = r, truth, lr_true = ns$.mbco_lr_at(prep, truth),
                    lr_a0 = ns$.mbco_branch_lr(prep)[["a0"]],
                    lr_b0 = ns$.mbco_branch_lr(prep)[["b0"]],
                    ci_lower = NA_real_, ci_upper = NA_real_)
  if (r <= width_reps) {
    fit <- tryCatch(pmed(fy, formula_m = fm, data = d, treatment = "X", mediator = "M",
                         method = "mbco"), error = function(e) NULL)
    if (!is.null(fit)) { out$ci_lower <- fit@ci_lower; out$ci_upper <- fit@ci_upper }
  }
  out
}
assign(".Random.seed", s, envir = globalenv())
res <- do.call(rbind, lapply(seq_len(reps), one))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, file.path(outdir, sprintf("null_task%03d.csv", task)), row.names = FALSE)
cat(sprintf("task %d (h_a=%g, h_b=%g, N=%d, cov=%s): %d reps, acceptance %.4f\n", task,
            pt$h_a, pt$h_b, pt$N, pt$cov, reps,
            mean(res$lr_true <= stats::qchisq(0.95, 1), na.rm = TRUE)))
