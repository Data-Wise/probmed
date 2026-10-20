# Read-only numerical and selected-stream verification. Run with one argument:
# oracle, i2, i6, or e. I2 runs from this checkout; I6/e run on Hopper.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L, args[[1L]] %in% c("oracle", "i2", "i6", "e"))
mode <- args[[1L]]

score <- function(a, b, vm, vy, delta = 1) {
  signal <- delta * a * b
  variance <- 2 * (b * b * vm + vy)
  if (any(!is.finite(c(signal, variance))) || variance < 0) return(NA_real_)
  if (variance == 0) return((sign(signal) + 1) / 2)
  stats::pnorm(signal / sqrt(variance))
}

oracle <- function(a, b, vm, vy, delta = 1) {
  if (vy == 0 && vm == 0) return((sign(a * b * delta) + 1) / 2)
  if (vy == 0 && b == 0) return(.5)
  if (vy == 0) return(stats::pnorm(sign(b) * a * delta / sqrt(2 * vm)))
  stats::integrate(function(z) {
    difference <- a * delta + sqrt(2 * vm) * z
    stats::pnorm(b * difference / sqrt(2 * vy)) * stats::dnorm(z)
  }, -Inf, Inf, rel.tol = 1e-12, abs.tol = 0)$value
}

exact_fit <- function(d) {
  if (!identical(names(d), c("X", "M", "Y")) ||
        !all(vapply(d, is.numeric, logical(1))) ||
        !all(is.finite(as.matrix(d))) || nrow(d) <= 3L ||
        !all(sort(unique(d$X)) == c(0, 1))) return(NA_real_)
  n <- nrow(d)
  mediator <- .lm.fit(cbind(1, d$X), d$M)
  outcome <- .lm.fit(cbind(1, d$X, d$M), d$Y)
  if (mediator$rank != 2L || outcome$rank != 3L) return(NA_real_)
  score(mediator$coefficients[[2L]], outcome$coefficients[[3L]],
        sum(mediator$residuals^2) / (n - 2L),
        sum(outcome$residuals^2) / (n - 3L))
}

guarded_interval <- function(point, boot) {
  failed <- sum(!is.finite(boot))
  ci <- if (is.finite(point) && failed == 0L)
    unname(stats::quantile(boot, c(.025, .975), names = FALSE))
  else c(NA_real_, NA_real_)
  c(lower = ci[[1L]], upper = ci[[2L]], failed = failed)
}

if (mode == "oracle") {
  fixtures <- rbind(
    c(.5, .4, 1, 1, 1), c(-.5, .4, 1, 2, 1),
    c(.5, -.4, 2, .5, -2), c(.5, .4, 2, .5, -2),
    c(0, .4, 1, 1, 1), c(.5, 0, 1, 1, 1),
    c(2.5, 1.5, 1, 2, 1), c(.5, .4, 0, 1, 1)
  )
  difference <- apply(fixtures, 1L, function(x) {
    abs(do.call(score, as.list(x)) - do.call(oracle, as.list(x)))
  })
  stopifnot(max(difference) <= 1e-8,
            score(1, 1, 0, 0) == 1, score(-1, 1, 0, 0) == 0,
            score(0, 1, 0, 0) == .5, is.na(score(1, 1, NA_real_, 1)))
  set.seed(20261020)
  n <- 300L
  x <- stats::rbinom(n, 1L, .5)
  group <- factor(sample(c("a", "b"), n, TRUE))
  m <- .5 * x + .7 * (group == "b") + stats::rnorm(n)
  y <- .3 * x + .4 * m - .2 * (group == "b") + stats::rnorm(n)
  fit <- function(outcome) {
    fm <- stats::lm(m ~ x + group)
    fy <- stats::lm(outcome ~ x + m + group)
    score(coef(fm)[["x"]], coef(fy)[["m"]], stats::sigma(fm)^2,
          stats::sigma(fy)^2)
  }
  stopifnot(abs(fit(y) - fit(y + 5 * x + 2 * (group == "b"))) < 1e-12)
  d <- data.frame(X = x, M = m, Y = y)
  stopifnot(is.finite(exact_fit(d)))
  bad <- d
  bad$M[[1L]] <- NA_real_
  stopifnot(is.na(exact_fit(bad)))
  bad <- d
  bad$X <- 0
  stopifnot(is.na(exact_fit(bad)))
  bad <- d
  bad$group <- group
  stopifnot(is.na(exact_fit(bad)))
  failed_interval <- guarded_interval(exact_fit(d), c(exact_fit(d), exact_fit(bad)))
  stopifnot(failed_interval[["failed"]] == 1L,
            is.na(failed_interval[["lower"]]),
            is.na(failed_interval[["upper"]]))
  simulated <- .5 * .4 + .4 * stats::rnorm(1000000L, sd = sqrt(2)) +
    stats::rnorm(1000000L, sd = sqrt(2))
  simulation_error <- abs(mean(simulated > 0) - score(.5, .4, 1, 1))
  stopifnot(simulation_error < .002)
  cat(sprintf("oracle max_abs_error=%.3g; simulation_abs_error=%.4g; guards=PASS\n",
              max(difference), simulation_error))
} else if (mode == "e") {
  stopifnot(length(args) == 3L)
  source(file.path(args[[2L]], "config.R"))
  cells <- expand.grid(N = c(50L, 100L, 250L, 500L, 1000L),
                       scenario = c("s1", "s6", "s7", "s8"),
                       stringsAsFactors = FALSE)
  errors <- numeric(nrow(cells))
  for (cell in seq_len(nrow(cells))) {
    task <- (cell - 1L) * 30L + 1L
    scenario <- cells$scenario[[cell]]
    n <- cells$N[[cell]]
    saved <- read.csv(file.path(args[[3L]], sprintf("e_%03d.csv", task)))
    saved <- saved[saved$rep %in% 1:2, ]
    stopifnot(nrow(saved) == 2L, identical(sort(saved$rep), 1:2),
              all(saved$boot_attempted == 2000L), all(saved$boot_failed == 0L))
    RNGkind("L'Ecuyer-CMRG")
    set.seed(20261021L + task)
    got <- parallel::mclapply(1:200, function(i) {
      if (i > 2L) return(NULL)
      d <- SCENARIOS[[scenario]]$gen(n)
      point <- exact_fit(d)
      boots <- replicate(2000L, exact_fit(d[sample.int(n, n, TRUE), , drop = FALSE]))
      c(rep = i, estimate = point, guarded_interval(point, boots))
    }, mc.cores = 8L)
    got <- do.call(rbind, got[1:2])
    got <- got[order(got[, "rep"]), , drop = FALSE]
    saved <- saved[order(saved$rep), , drop = FALSE]
    stopifnot(all(got[, "failed"] == 0L))
    saved_values <- as.matrix(saved[, c("estimate", "ci_lower", "ci_upper")])
    errors[[cell]] <- max(abs(got[, c("estimate", "lower", "upper")] - saved_values))
    if (errors[[cell]] > 1e-12)
      stop(sprintf("e %s N=%d task=%d error=%.17g",
                   scenario, n, task, errors[[cell]]))
  }
  cat(sprintf("e selected replay: %d cells x 2 reps; max_abs_error=%.3g\n",
              nrow(cells), max(errors)))
} else {
  stopifnot(length(args) == 3L)
  source(file.path(args[[2L]], "config.R"))
  source(file.path(args[[2L]], "simulation_functions.R"))
  source(file.path(args[[2L]], "mbco_pmed.R"))
  cells <- expand.grid(N = SCALES$cluster$N, scenario = names(SCENARIOS),
                       stringsAsFactors = FALSE)
  scenarios <- if (mode == "i2") {
    c("s1", "s2", "s6", "s7", "s8")
  } else {
    c("s1", "s6", "s7", "s8")
  }
  i2 <- if (mode == "i2") read.csv(args[[3L]]) else NULL
  if (mode == "i2") {
    keys <- paste(i2$scenario, i2$N, i2$rep)
    grid <- expand.grid(rep = seq_len(2000L), N = SCALES$cluster$N,
                        scenario = scenarios, stringsAsFactors = FALSE)
    expected <- paste(grid$scenario, grid$N, grid$rep)
    stopifnot(identical(sort(keys), sort(expected)),
              all(i2$boot_attempted == 1000L), all(i2$boot_failed == 0L),
              all(is.finite(as.matrix(i2[, c("exact_estimate", "ci_lower", "ci_upper")]))),
              all(i2$ci_lower <= i2$ci_upper))
    cat(sprintf("I2 grid: %d unique expected rows; finite intervals; failures=0\n",
                nrow(i2)))
  } else {
    active <- cells[cells$scenario %in% c("s1", "s3", "s4", "s6", "s7", "s8"), ]
    rows <- 0L
    for (cell_idx in which(cells$scenario %in% active$scenario)) {
      for (chunk in seq_len(10L)) {
        task_id <- (cell_idx - 1L) * 10L + chunk
        path <- file.path(args[[3L]], sprintf("i6_%d.csv", task_id))
        stopifnot(file.exists(path))
        part <- read.csv(path)
        expected_rep <- ((chunk - 1L) * 200L + 1L):(chunk * 200L)
        stopifnot(nrow(part) == 200L,
                  identical(sort(part$rep), expected_rep),
                  all(part$scenario == cells$scenario[[cell_idx]]),
                  all(part$N == cells$N[[cell_idx]]),
                  all(part$boot_attempted == 2000L))
        if (cells$scenario[[cell_idx]] %in% scenarios) {
          fields <- c("exact_estimate", "ci_lower_e", "ci_upper_e",
                      "ci_lower_d", "ci_upper_d")
          stopifnot(all(is.finite(as.matrix(part[, fields]))),
                    all(part$ci_lower_e <= part$ci_upper_e),
                    all(part$boot_failed_d == 0L))
        } else {
          stopifnot(all(is.finite(as.matrix(part[, c("ci_lower_b", "ci_upper_b")]))),
                    all(part$ci_lower_b <= part$ci_upper_b),
                    all(part$boot_failed_b == 0L))
        }
        rows <- rows + nrow(part)
      }
    }
    stopifnot(rows == 60000L)
    cat(sprintf("I6 grid: %d expected rows across 300 parts; finite active intervals; failures=0\n",
                rows))
  }
  differences <- numeric()
  for (scenario in scenarios) for (n in SCALES$cluster$N) {
    cell <- which(cells$scenario == scenario & cells$N == n)
    task <- (cell - 1L) * 10L + 1L
    saved <- if (mode == "i2") {
      i2[i2$scenario == scenario & i2$N == n & i2$rep %in% 1:2, ]
    } else {
      path <- file.path(args[[3L]], sprintf("i6_%d.csv", task))
      subset(read.csv(path), rep %in% 1:2)
    }
    stopifnot(nrow(saved) == 2L, identical(sort(saved$rep), 1:2))
    boot_n <- if (mode == "i2") 1000L else 2000L
    seed_base <- if (mode == "i2") 20261016L else 20261017L
    replay_one <- function(i) {
      d <- SCENARIOS[[scenario]]$gen(n)
      invisible(estimate_all(d, "gaussian", SCALES$cluster$K))
      u <- crn_uniforms(50000L)
      invisible(pmed_point(d, "gaussian", 50000L, u))
      invisible(pmed_point(d, "gaussian", 100L * SCALES$cluster$K))
      idx <- replicate(boot_n, sample.int(n, n, replace = TRUE))
      boots <- apply(idx, 2L, function(ix) exact_fit(d[ix, , drop = FALSE]))
      point <- exact_fit(d)
      c(rep = i, estimate = point, guarded_interval(point, boots))
    }
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed_base + task)
    got <- parallel::mclapply(1:200, function(i) {
      if (i > 2L) return(NULL)
      replay_one(i)
    }, mc.cores = 8L)
    got <- do.call(rbind, got[1:2])
    got <- got[order(got[, "rep"]), , drop = FALSE]
    saved <- saved[order(saved$rep), , drop = FALSE]
    target <- if (mode == "i2")
      cbind(saved$exact_estimate, saved$ci_lower, saved$ci_upper)
    else cbind(saved$exact_estimate, saved$ci_lower_e, saved$ci_upper_e)
    stopifnot(all(got[, "failed"] == 0L), all(is.finite(target)))
    err <- max(abs(got[, c("estimate", "lower", "upper")] - target))
    if (err > 1e-12) {
      stop(sprintf("%s %s N=%d task=%d error=%.17g",
                   mode, scenario, n, task, err))
    }
    differences <- c(differences, err)
  }
  cat(sprintf("%s selected replay: %d cells x 2 reps; max_abs_error=%.3g\n",
              mode, length(differences), max(differences)))
}
