# Independent, read-only audit of the frozen exact-Gaussian Hopper output.
# Usage: Rscript --vanilla audit_exact_gaussian_hopper.R /absolute/run/path

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, startsWith(args[[1L]], "/"))
root <- args[[1L]]
manifest <- read.csv(file.path(root, "code", "manifest.csv"))
saved <- read.csv(file.path(root, "summary.csv"))

source_files <- c(
  runner_sha256 = "runner.R",
  kernel_sha256 = "kernel.R",
  config_sha256 = "config.R",
  analyze_sha256 = "analyze.R",
  submit_sha256 = "submit_exact.sh"
)
sha256 <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE)
  stopifnot(length(output) == 1L, grepl("^[0-9a-f]{64}  ", output))
  substr(output, 1L, 64L)
}
stopifnot(nrow(manifest) == 600L, identical(manifest$task, 1:600),
          all(manifest$boot_n == 2000L),
          all(manifest$seed_base == 20261021L))
for (column in names(source_files)) {
  actual <- sha256(file.path(root, "code", source_files[[column]]))
  stopifnot(all(manifest[[column]] == actual))
}

cells <- expand.grid(
  N = c(50L, 100L, 250L, 500L, 1000L),
  scenario = c("s1", "s6", "s7", "s8"),
  stringsAsFactors = FALSE
)
expected_cell <- (manifest$task - 1L) %/% 30L + 1L
expected_chunk <- (manifest$task - 1L) %% 30L + 1L
stopifnot(identical(manifest$scenario, cells$scenario[expected_cell]),
          identical(manifest$N, cells$N[expected_cell]),
          identical(manifest$rep_first, (expected_chunk - 1L) * 200L + 1L),
          identical(manifest$rep_last, expected_chunk * 200L))

parts_dir <- file.path(root, "parts")
expected_files <- file.path(parts_dir, sprintf("e_%03d.csv", 1:600))
actual_files <- list.files(parts_dir, pattern = "^e_[0-9]{3}\\.csv$", full.names = TRUE)
actual_sidecars <- list.files(parts_dir, pattern = "^e_[0-9]{3}\\.csv\\.sha256$", full.names = TRUE)
stopifnot(setequal(actual_files, expected_files),
          setequal(actual_sidecars, paste0(expected_files, ".sha256")))
hash_lines <- readLines(file.path(root, "part_hashes.sha256"), warn = FALSE)
hashes <- substr(hash_lines, 1L, 64L)
paths <- substring(hash_lines, 66L)
stopifnot(length(hash_lines) == 600L,
          all(grepl("^[0-9a-f]{64} /", hash_lines)),
          identical(paths, expected_files))
for (i in seq_along(expected_files)) {
  stopifnot(identical(readLines(paste0(expected_files[[i]], ".sha256"),
                            warn = FALSE), hashes[[i]]))
}

expected_columns <- c(
  "task", "scenario", "N", "rep", "estimate", "ci_lower", "ci_upper",
  "boot_attempted", "boot_failed", "seed_base", "runner_sha256",
  "kernel_sha256", "config_sha256"
)
parts <- lapply(seq_along(expected_files), function(i) {
  x <- read.csv(expected_files[[i]])
  m <- manifest[i, ]
  stopifnot(identical(names(x), expected_columns), nrow(x) == 200L,
            all(x$task == m$task), all(x$scenario == m$scenario),
            all(x$N == m$N), identical(x$rep, m$rep_first:m$rep_last),
            all(x$boot_attempted == 2000L),
            all(x$seed_base == 20261021L),
            all(x$runner_sha256 == m$runner_sha256),
            all(x$kernel_sha256 == m$kernel_sha256),
            all(x$config_sha256 == m$config_sha256),
            all(!is.na(x$boot_failed)),
            all(x$boot_failed >= 0L & x$boot_failed <= 2000L))
  complete <- is.finite(x$estimate) & is.finite(x$ci_lower) &
    is.finite(x$ci_upper) & x$boot_failed == 0L
  stopifnot(all(x$estimate[is.finite(x$estimate)] >= 0 &
                  x$estimate[is.finite(x$estimate)] <= 1),
            all(x$ci_lower[complete] >= 0 & x$ci_lower[complete] <=
                  x$ci_upper[complete] & x$ci_upper[complete] <= 1),
            all(is.na(x$ci_lower[!complete]) & is.na(x$ci_upper[!complete])))
  x
})
x <- do.call(rbind, parts)
keys <- paste(x$scenario, x$N, x$rep, sep = ":")
stopifnot(nrow(x) == 120000L, !anyDuplicated(keys),
          nrow(saved) == 20L,
          !anyDuplicated(paste(saved$scenario, saved$N, sep = ":")))

rows <- lapply(seq_len(nrow(cells)), function(i) {
  cell <- cells[i, ]
  z <- x[x$scenario == cell$scenario & x$N == cell$N, ]
  s <- saved[saved$scenario == cell$scenario & saved$N == cell$N, ]
  stopifnot(nrow(z) == 6000L, identical(z$rep, 1:6000), nrow(s) == 1L)
  truth <- if (cell$scenario %in% c("s6", "s7")) 0.5 else s$truth
  complete <- is.finite(z$estimate) & is.finite(z$ci_lower) &
    is.finite(z$ci_upper) & z$boot_failed == 0L
  covered <- complete & z$ci_lower <= truth & truth <= z$ci_upper
  covered[is.na(covered)] <- FALSE
  k <- sum(covered)
  lower <- if (k == 0L) 0 else qbeta(.05, k, 6001L - k)
  upper <- if (k == 6000L) 1 else qbeta(.95, k + 1L, 6000L - k)
  lower95 <- if (k == 0L) 0 else qbeta(.025, k, 6001L - k)
  upper95 <- if (k == 6000L) 1 else qbeta(.975, k + 1L, 6000L - k)
  width <- if (any(complete)) {
    mean(z$ci_upper[complete] - z$ci_lower[complete])
  } else NA_real_
  mean_estimate <- if (any(is.finite(z$estimate))) {
    mean(z$estimate[is.finite(z$estimate)])
  } else NA_real_
  checks <- c(
    attempted = nrow(z), incomplete = sum(!complete),
    failed_boot_fits = sum(z$boot_failed), covered = k,
    coverage = k / 6000, coverage_95_lo = lower95,
    coverage_95_hi = upper95, coverage_90_lo = lower,
    coverage_90_hi = upper, mean_width_complete = width,
    mean_estimate_finite = mean_estimate
  )
  for (field in names(checks)) {
    stopifnot(isTRUE(all.equal(unname(checks[[field]]), s[[field]],
                               tolerance = 1e-10)))
  }
  data.frame(scenario = cell$scenario, N = cell$N, as.list(checks))
})
audited <- do.call(rbind, rows)
gate <- audited$scenario %in% c("s6", "s7") & audited$N >= 250L
gate_rows <- audited[gate, ]
gate_rows$mcse <- sqrt(gate_rows$coverage * (1 - gate_rows$coverage) /
                         gate_rows$attempted)
gate_rows$pass <- gate_rows$incomplete == 0L &
  gate_rows$coverage_90_lo > .93 & gate_rows$coverage_90_hi < .97
stopifnot(nrow(gate_rows) == 6L)
print(gate_rows, row.names = FALSE, digits = 8)
cat("all_six_pass:", all(gate_rows$pass), "\n")
cat("rows:", nrow(x), "chunks:", length(parts),
    "unique_keys:", length(unique(keys)), "\n")
cat("all_cells_incomplete:", sum(audited$incomplete),
    "all_cells_failed_boot_fits:", sum(audited$failed_boot_fits), "\n")
cat("s6_N50_coverage:", audited$coverage[audited$scenario == "s6" &
                                      audited$N == 50L], "\n")
cat("summary_sha256:", sha256(file.path(root, "summary.csv")), "\n")
