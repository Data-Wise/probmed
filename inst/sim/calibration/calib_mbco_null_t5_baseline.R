# Phase 8 T5 baseline and branch-order audit using frozen T2 replications.

reconstruct_t5_branches <- function(d) {
  df_m <- d$N - 2L - as.integer(d$cov)
  df_y <- d$N - 3L - as.integer(d$cov)
  t_a <- sqrt(df_m * expm1(d$lr_a0 / d$N))
  t_b <- sqrt(df_y * expm1(d$lr_b0 / d$N))
  a_first <- t_a <= t_b
  data.frame(
    t_a = t_a, t_b = t_b, df_m = df_m, df_y = df_y, a_first = a_first,
    branch_df = ifelse(a_first, df_m, df_y),
    lr_disagrees = a_first != (d$lr_a0 <= d$lr_b0)
  )
}

t5_disposition <- function(cells) {
  if ("FAIL" %in% cells) return("FAIL")
  if (all(cells == "PASS")) return("PASS")
  "INCONCLUSIVE"
}

t5_cell_gates <- function(observed) {
  observed$accepted <- as.integer(round(observed$reps * observed$coverage))
  observed$minimum_gate <- ifelse(
    observed$ci_lower >= 0.94, "PASS",
    ifelse(observed$ci_upper < 0.94, "FAIL", "INCONCLUSIVE")
  )
  origin <- observed$h_a == 0 & observed$h_b == 0
  observed$origin_gate <- NA_character_
  observed$origin_gate[origin] <- ifelse(
    observed$ci_upper[origin] <= 0.97, "PASS",
    ifelse(observed$ci_lower[origin] > 0.97, "FAIL", "INCONCLUSIVE")
  )
  observed
}

summarize_t5_gates <- function(observed) {
  do.call(rbind, lapply(split(observed, observed$rule), function(cell) {
    minimum <- t5_disposition(cell$minimum_gate)
    origin <- t5_disposition(cell$origin_gate[!is.na(cell$origin_gate)])
    data.frame(rule = cell$rule[1], minimum, origin,
               overall = t5_disposition(c(minimum, origin)))
  }))

}

audit_t5_baseline <- function(d, baseline, score_rules) {
  rep_key <- c("h_a", "h_b", "N", "cov", "rep")
  score_key <- c("rule", "h_a", "h_b", "N", "cov")
  effects <- c(0, 0.5, 1, 2, 4)
  lrs <- as.matrix(d[c("lr_true", "lr_a0", "lr_b0")])
  stopifnot(
    nrow(d) == 200000L,
    is.logical(d$cov), !anyNA(d$cov), !anyDuplicated(d[rep_key]),
    all(d$h_a %in% effects), all(d$h_b %in% effects),
    all(d$N %in% c(100L, 1000L)), all(d$rep %in% seq_len(2000L)),
    all(is.finite(lrs) & lrs >= 0)
  )

  observed <- score_rules(d)
  comparison <- merge(
    observed, baseline, by = score_key, suffixes = c(".observed", ".saved")
  )
  stopifnot(
    nrow(observed) == 400L, nrow(baseline) == 400L,
    !anyDuplicated(observed[score_key]), !anyDuplicated(baseline[score_key]),
    nrow(comparison) == 400L,
    all(comparison$reps.observed == 2000L),
    all(comparison$reps.observed == comparison$reps.saved),
    all(round(comparison$coverage.observed * comparison$reps.observed) ==
          round(comparison$coverage.saved * comparison$reps.saved)),
    all(abs(comparison$coverage.observed - comparison$coverage.saved) <= 1e-12),
    all(abs(comparison$ci_lower.observed - comparison$ci_lower.saved) <= 1e-12),
    all(abs(comparison$ci_upper.observed - comparison$ci_upper.saved) <= 1e-12)
  )

  observed <- t5_cell_gates(observed)

  list(
    scores = observed, gates = summarize_t5_gates(observed),
    branch_disagreements = sum(reconstruct_t5_branches(d)$lr_disagrees)
  )
}

write_t5_scores <- function(scores, output_path) {
  staged <- tempfile(".t5-scores-", tmpdir = dirname(output_path))
  on.exit(unlink(staged), add = TRUE)
  utils::write.csv(scores, staged, row.names = FALSE)
  # Replace the directory entry so existing hard links keep their original data.
  if (!file.rename(staged, output_path)) {
    stop("Could not replace T5 score output: ", output_path)
  }
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3L) {
    stop("Usage: Rscript calib_mbco_null_t5_baseline.R replications.csv.gz scores.csv output.csv")
  }
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  source(file.path(dirname(script), "calib_mbco_null_score.R"))
  paths <- normalizePath(args, mustWork = FALSE)
  if (paths[3L] %in% paths[1:2]) stop("Use a separate output path for the baseline audit.")
  result <- audit_t5_baseline(utils::read.csv(args[1L]), utils::read.csv(args[2L]), score_rules)
  write_t5_scores(result$scores, args[3L])
  print(result$gates, row.names = FALSE)
  cat("LR vs |t| branch-order disagreements:", result$branch_disagreements, "\n")
}
