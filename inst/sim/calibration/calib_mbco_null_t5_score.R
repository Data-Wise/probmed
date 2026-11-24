# Offline T5 audit of the frozen MBCO lattice, including the four T4 baselines.
# Source the T4 scorer and T5 baseline helpers before calling audit_t5_prototypes().

t5_checkerboard_reject <- function(u_a, u_b, alpha) {
  (u_a >= 1 - alpha & u_b >= 1 - alpha) |
    (floor(u_a / alpha) == floor(u_b / alpha))
}

t5_prototype_acceptance <- function(d) {
  alpha <- 0.05
  c_hi <- stats::qchisq(1 - alpha, 1)
  c_lo <- stats::qchisq(1 - sqrt(alpha), 1)
  z_hi <- sqrt(c_hi)
  null <- d$h_a == 0 | d$h_b == 0
  chi2_accept <- d$lr_true <= c_hi
  null_only <- function(accept) ifelse(null, accept, chi2_accept)

  branches <- reconstruct_t5_branches(d)
  u_a_t <- 2 * stats::pt(branches$t_a, df = branches$df_m) - 1
  u_b_t <- 2 * stats::pt(branches$t_b, df = branches$df_y) - 1
  z_a <- sqrt(d$lr_a0)
  z_b <- sqrt(d$lr_b0)
  z_max <- pmax(z_a, z_b)
  u_a_lr <- 2 * stats::pnorm(z_a) - 1
  u_b_lr <- 2 * stats::pnorm(z_b) - 1

  # Historical all-p* LR heuristic. No t-to-LR cutoff mapping is scored.
  u_max_lr <- pmax(u_a_lr, u_b_lr)
  bin_edge <- floor(u_max_lr / alpha) * alpha
  mc21_lr_cut <- ifelse(u_max_lr < 1 - alpha,
                        stats::qnorm((1 + bin_edge) / 2)^2, c_hi)

  # The prototype adds a strip with origin probability alpha - alpha^2.
  strip_start <- 0.25
  strip_probability <- function(width) {
    stats::integrate(function(z) {
      8 * stats::dnorm(z) * (stats::pnorm(z) - stats::pnorm(pmax(0, z - width)))
    }, lower = strip_start, upper = z_hi)$value
  }
  strip_width <- stats::uniroot(
    function(width) strip_probability(width) - (alpha - alpha^2),
    interval = c(0.001, 1), tol = 1e-9
  )$root
  strip_cut <- ifelse(z_max >= strip_start & z_max < z_hi,
                      pmax(0, z_max - strip_width)^2, c_hi)

  # Analytical JS analogues. Interpolate the squared cutoff on the original grid.
  shifts <- seq(0, 8, by = 0.01)
  js_cutoffs <- vapply(shifts, function(mu) {
    if (mu == 0) return(sqrt(c_lo))
    tail_error <- function(cutoff) {
      2 * stats::pnorm(-cutoff) *
        (stats::pnorm(-cutoff - mu) + stats::pnorm(-cutoff + mu)) - alpha
    }
    stats::uniroot(tail_error, c(sqrt(c_lo), z_hi + 0.01), tol = 1e-10)$root
  }, numeric(1))^2
  js_cut <- stats::approxfun(shifts, js_cutoffs, rule = 2)
  logn_accept <- d$lr_true <= js_cut(ifelse(z_max > sqrt(2 * log(d$N)), z_max, 0))
  logn2_accept <- d$lr_true <= js_cut(ifelse(z_max > 2 * sqrt(log(d$N)), z_max, 0))
  plugin_accept <- d$lr_true <= js_cut(z_max)
  strip_accept <- d$lr_true <= strip_cut

  list(
    MC21_t_05 = null_only(!t5_checkerboard_reject(u_a_t, u_b_t, alpha)),
    MC21_LR_05 = null_only(!t5_checkerboard_reject(u_a_lr, u_b_lr, alpha)),
    MC21_LR_all = d$lr_true <= mc21_lr_cut,
    vGvG_strip_05 = null_only(strip_accept),
    vGvG_strip_all = strip_accept,
    HSX24_logN_05 = null_only(logn_accept),
    HSX24_logN_all = logn_accept,
    HSX24_2logN_05 = null_only(logn2_accept),
    HSX24_2logN_all = logn2_accept,
    HSX24_plugin_05 = null_only(plugin_accept),
    HSX24_plugin_all = plugin_accept
  )
}

audit_t5_prototypes <- function(d, baseline, score_rules) {
  audit <- audit_t5_baseline(d, baseline, score_rules)
  indices <- split(seq_len(nrow(d)), list(d$h_a, d$h_b, d$N, d$cov), drop = TRUE)
  keys <- d[vapply(indices, function(i) i[1L], integer(1)), c("h_a", "h_b", "N", "cov")]
  accept <- t5_prototype_acceptance(d)
  stopifnot(all(lengths(indices) == 2000L), !anyNA(unlist(accept, use.names = FALSE)))
  prototypes <- lapply(names(accept), function(rule) {
    accepted <- vapply(indices, function(i) sum(accept[[rule]][i]), integer(1))
    ci <- vapply(accepted, function(n) stats::binom.test(n, 2000L)$conf.int, numeric(2))
    data.frame(rule, keys, reps = 2000L, coverage = accepted / 2000L,
               ci_lower = ci[1L, ], ci_upper = ci[2L, ], row.names = NULL)
  })
  scores <- rbind(audit$scores, t5_cell_gates(do.call(rbind, prototypes)))
  scores <- scores[order(scores$rule, scores$N, scores$cov, scores$h_a, scores$h_b,
                         method = "radix"), ]
  rownames(scores) <- NULL
  stopifnot(nrow(scores) == 1500L,
            !anyDuplicated(scores[c("rule", "h_a", "h_b", "N", "cov")]))
  audit$scores <- scores
  audit$gates <- summarize_t5_gates(scores)
  audit
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) {
    stop("Usage: Rscript calib_mbco_null_t5_score.R replications.csv.gz output.csv")
  }
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  source(file.path(dirname(script), "calib_mbco_null_score.R"))
  source(file.path(dirname(script), "calib_mbco_null_t5_baseline.R"))
  baseline <- file.path(dirname(args[1L]), "scores.csv")
  input_paths <- normalizePath(c(args[1L], baseline), mustWork = TRUE)
  output_path <- file.path(normalizePath(dirname(args[2L]), mustWork = TRUE),
                           basename(args[2L]))
  if (normalizePath(output_path, mustWork = FALSE) %in% input_paths) {
    stop("Use a separate output path for the prototype audit.")
  }
  result <- audit_t5_prototypes(utils::read.csv(args[1L]), utils::read.csv(baseline), score_rules)
  write_t5_scores(result$scores, output_path)
  print(result$gates, row.names = FALSE)
  cat("LR vs |t| branch-order disagreements:", result$branch_disagreements, "\n")
}
