#!/usr/bin/env Rscript
# G2c: audit a proposed fresh seed block against every earlier seed range used in this project.
# Usage: Rscript scripts/audit-seed-ranges.R [pr55_scripts_dir]
#   pr55_scripts_dir defaults to /private/tmp/probmed-dev-main-review/scripts (the PR #55 checkout).
# Proposed block: seed = 1,000,000,000 + cell * 10000 + replicate, cell 1..40, replicate 1..5000.
# NOT registered or frozen here; the block is frozen only by T5 after the author's gates exist.
args <- commandArgs(TRUE)
pr55 <- if (length(args)) args[1] else "/private/tmp/probmed-dev-main-review/scripts"
prop <- c(lo = 1e9 + 1 * 10000 + 1, hi = 1e9 + 40 * 10000 + 5000)

# Earlier ranges, from the formulas read in the repo / PR #55 scripts (source in the last column).
# Widths are upper bounds: cell <= 40 and replicate <= 5000 for the cell*10000+r studies; hopper i4/i6 use
# base + task*1100000 + i*5000 + j*2 + k with task <= 400 assumed (arrays used are 1-300).
earlier <- data.frame(
  name = c("PR55 observed (20261020 + cell*1e4 + r)", "PR55 bias diagnostics (20261020 + 20000 + r)",
    "PR55 boundary-correction (20361020 + ...)", "PR55 joint-contrast (20461020 + ...)",
    "PR55 selection (20561020 + ...)", "hopper i4 (20261016 + task*1.1e6 + ...)",
    "hopper i6 (20261017 + task*1.1e6 + ...)", "calibration (20261014-20261015)",
    "exact Gaussian (20261019-20261021 + ...)", "inst/sim seed 1000000 family", "pmedw 424242", "test 20281186",
    "this branch scripts (20261022, 20261023 + n, small ints)"),
  lo = c(20261020 + 10000, 20261020 + 20000, 20361020 + 10000, 20461020 + 10000, 20561020 + 10000,
    20261016, 20261017, 20261014, 20261019, 1e6, 424242, 20281186, 1),
  hi = c(20261020 + 70000 + 5000, 20261020 + 20000 + 5000, 20361020 + 70000 + 5000, 20461020 + 70000 + 5000,
    20561020 + 70000 + 5000, 20261016 + 400 * 1100000 + 1100000, 20261017 + 400 * 1100000 + 1100000,
    20261015 + 1e6, 20261021 + 1e6, 1e7, 424242 + 1e6, 20281186, 20261023 + 1000))
overlap <- function(lo, hi) earlier$name[earlier$lo <= hi & lo <= earlier$hi]

# Controls: a planted overlapping block must be flagged; a planted disjoint block must not.
stopifnot(length(overlap(20461020 + 10001, 20461020 + 10005)) > 0, length(overlap(5e8, 5e8 + 10)) == 0)
hit <- overlap(prop[["lo"]], prop[["hi"]])
cat(sprintf("Proposed block [%.0f, %.0f]; max 32-bit integer seed %.0f: %s\n", prop[["lo"]], prop[["hi"]], 2^31 - 1,
  if (prop[["hi"]] < 2^31 - 1) "fits" else "DOES NOT FIT"))
cat(sprintf("Overlap with %d earlier ranges: %s\n", nrow(earlier), if (length(hit)) paste(hit, collapse = "; ") else "none"))

# Independent discovery: every integer literal >= 6 digits on a line mentioning seed, in this repo and PR #55 scripts.
# The PR #55 scripts are required evidence: fail loudly if the directory or the expected scripts are absent
# (Codex review 2026-10-21, finding 1), rather than skipping them and printing AUDIT PASS.
expected <- c("validate-wasserstein-observed.R", "validate-wasserstein-boundary-correction.R",
  "validate-wasserstein-joint-contrast.R", "validate-wasserstein-boundary-selection.R")
if (!dir.exists(pr55) || !all(expected %in% list.files(pr55))) {
  cat(sprintf("AUDIT FAIL: PR #55 scripts missing in '%s' (need %s)\n", pr55, paste(expected, collapse = ", ")))
  quit(status = 2)
}
files <- c(list.files(c("scripts", "inst", "R", "tests"), pattern = "[.](R|sh|dcf)$", recursive = TRUE, full.names = TRUE),
  list.files(pr55, pattern = "[.]R$", full.names = TRUE))
found <- list()
for (f in files) {
  ln <- readLines(f, warn = FALSE)
  for (i in grep("seed|SEED", ln)) {
    m <- regmatches(ln[i], gregexpr("[0-9]{6,10}", ln[i]))[[1]]
    for (v in as.numeric(m)) found[[length(found) + 1L]] <- data.frame(file = f, line = i, value = v)
  }
}
found <- do.call(rbind, found)
inside <- found[found$value >= prop[["lo"]] & found$value <= prop[["hi"]], ]
cat(sprintf("Discovery: %d seed-line integer literals in %d files; inside the proposed block: %d\n", nrow(found), length(files), nrow(inside)))
if (nrow(inside)) print(inside, row.names = FALSE)
cat(if (length(hit) == 0 && nrow(inside) == 0 && prop[["hi"]] < 2^31 - 1) "AUDIT PASS\n" else "AUDIT FAIL\n")
cat("Limit: only seed literals in this repo and the PR #55 scripts directory are covered; hopper jobs run from other\n",
  "checkouts (for example ~/gauge_boot_postfix/) are not, and are blocked here only by the stated 1e9 distance.\n", sep = "")
