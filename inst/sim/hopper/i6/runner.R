args <- commandArgs(trailingOnly = TRUE)
task <- as.integer(args[[1]])
out_file <- args[[2]]
src <- Sys.getenv("I6_PMED_SRC", ".")
source(file.path(src, "config.R"))
source(file.path(src, "simulation_functions.R"))
source(file.path(src, "mbco_pmed.R"))

sc <- SCALES[["cluster"]]
cells <- expand.grid(N = sc$N, scenario = names(SCENARIOS), stringsAsFactors = FALSE)
N_CHUNKS <- 10L
stopifnot(task >= 1L, task <= nrow(cells) * N_CHUNKS)
cell_idx <- (task - 1L) %/% N_CHUNKS + 1L
chunk <- (task - 1L) %% N_CHUNKS + 1L
sn <- cells$scenario[cell_idx]; N <- cells$N[cell_idx]
per <- sc$reps %/% N_CHUNKS
reps <- ((chunk - 1L) * per + 1L):(chunk * per)
pilot_n <- as.integer(Sys.getenv("I6_PILOT_REPS", "0"))
if (pilot_n > 0L) reps <- head(reps, pilot_n)
boot_n <- as.integer(Sys.getenv("I6_BOOT_N", as.character(sc$n_boot)))
K_POINT <- sc$K
K_BOOT <- as.integer(Sys.getenv("I6_K_BOOT", "50000"))
K_LARGE <- as.integer(Sys.getenv("I6_K_LARGE", if (SCENARIOS[[sn]]$family_y == "gaussian") "200000" else "1000000"))
base <- as.integer(Sys.getenv("I6_SEED_BASE", "20261017"))
N_CHECK <- as.integer(Sys.getenv("I6_N_CHECK", "20"))
cores <- as.integer(Sys.getenv("I6_CORES", "8"))
preflight <- identical(Sys.getenv("I6_PREFLIGHT", "0"), "1")

exact_pmed <- function(d) {
  n <- nrow(d); fm <- .lm.fit(cbind(1, d$X), d$M); fy <- .lm.fit(cbind(1, d$X, d$M), d$Y)
  a <- fm$coefficients[2L]; b <- fy$coefficients[3L]
  sm2 <- sum(fm$residuals^2) / (n - 2L); sy2 <- sum(fy$residuals^2) / (n - 3L)
  as.numeric(pnorm(a * b / sqrt(2 * (b^2 * sm2 + sy2))))
}
pmed_scores <- function(data, family_y, K) {
  fm <- lm(M ~ X, data)
  b0m <- coef(fm)["(Intercept)"]; am <- coef(fm)["X"]; sm <- sigma(fm)
  M0 <- rnorm(K, b0m, sm); M1 <- rnorm(K, b0m + am, sm)
  if (family_y == "gaussian") {
    fy <- lm(Y ~ X + M, data); cy <- coef(fy); sy <- sigma(fy)
    lin <- function(m) cy["(Intercept)"] + cy["X"] + cy["M"] * m
    Y0 <- lin(M0) + rnorm(K, 0, sy); Y1 <- lin(M1) + rnorm(K, 0, sy)
    return(as.numeric(Y1 > Y0))
  }
  fy <- suppressWarnings(glm(Y ~ X + M, data, family = binomial))
  cy <- coef(fy)
  p <- function(m) plogis(cy["(Intercept)"] + cy["X"] + cy["M"] * m)
  Y0 <- rbinom(K, 1, p(M0)); Y1 <- rbinom(K, 1, p(M1))
  as.numeric(Y1 > Y0) + 0.5 * as.numeric(Y1 == Y0)
}
score_interval <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  unname(quantile(x, c(.025, .975), names = FALSE))
}
one <- function(i) {
  scenario <- SCENARIOS[[sn]]; fy <- scenario$family_y
  d <- scenario$gen(N)
  invisible(estimate_all(d, fy, K_POINT))
  u <- crn_uniforms(K_BOOT)
  anchor <- as.numeric(pmed_point(d, fy, K_BOOT, u))
  point <- as.numeric(pmed_point(d, fy, 100L * K_POINT))
  idx <- replicate(boot_n, sample.int(N, N, replace = TRUE))
  exact <- if (fy == "gaussian") apply(idx, 2L, function(ix) exact_pmed(d[ix, , drop = FALSE])) else rep(NA_real_, boot_n)
  worker_kind <- RNGkind()
  worker_seed <- .Random.seed
  raw <- cv <- large <- crn <- rep(NA_real_, boot_n)
  var_mc <- rep(NA_real_, boot_n)
  check_z <- rep(NA_real_, boot_n)
  fail_a <- fail_b <- fail_d <- 0L
  if (!preflight) {
    for (j in seq_len(boot_n)) {
      db <- d[idx[, j], , drop = FALSE]
      old_kind <- RNGkind()
      set.seed(base + task * 1100000L + i * 5000L + j * 2L + 1L)
      scores <- tryCatch(pmed_scores(db, fy, K_POINT), error = function(e) NULL)
      if (fy == "gaussian") {
        # I6: Gaussian cells run d (fresh K_LARGE draws) against in-run exact (e)
        # Given the fit, Y1 > Y0 is iid Bernoulli(p) with p = exact_pmed(resample), so the K_LARGE-draw
        # mean is exactly Binomial(K_LARGE, p) / K_LARGE; the first N_CHECK resamples also run the real
        # pmed_scores() path and record its standardized deviation from p (should look N(0, 1)).
        set.seed(base + task * 1100000L + i * 5000L + j * 2L + 2L)
        large[j] <- stats::rbinom(1L, K_LARGE, exact[j]) / K_LARGE
        if (j <= N_CHECK) {
          set.seed(base + task * 1100000L + i * 5000L + j * 2L + 3L)
          chk <- tryCatch(pmed_scores(db, fy, K_LARGE), error = function(e) NULL)
          if (!is.null(chk) && length(chk)) check_z[j] <- (mean(chk) - exact[j]) / sqrt(exact[j] * (1 - exact[j]) / K_LARGE)
        }
      } else {
        # I6: binary cells are report-only and run b (CRN + control variate)
        z <- tryCatch(as.numeric(pmed_point(db, fy, K_BOOT, u)), error = function(e) NA_real_)
        if (is.finite(z)) crn[j] <- min(max(z - anchor + point, 0), 1) else fail_b <- fail_b + 1L
      }
      RNGkind(old_kind[1L], old_kind[2L], old_kind[3L])
    }
  }
  out <- list(
    scenario = sn, N = N, rep = i,
    ci_lower_a = NA_real_, ci_upper_a = NA_real_,
    ci_lower_b = NA_real_, ci_upper_b = NA_real_,
    ci_lower_c = NA_real_, ci_upper_c = NA_real_,
    ci_lower_d = NA_real_, ci_upper_d = NA_real_,
    ci_lower_e = NA_real_, ci_upper_e = NA_real_,
    boot_failed_a = fail_a, boot_failed_b = fail_b, boot_failed_d = fail_d,
    boot_attempted = boot_n, correction_factor = NA_real_,
    exact_estimate = if (fy == "gaussian") exact_pmed(d) else NA_real_
  )
  if (fy == "gaussian") {
    q <- score_interval(exact); out$ci_lower_e <- q[1]; out$ci_upper_e <- q[2]
  }
  if (!preflight) {
    if (fy == "gaussian") {
      q <- score_interval(large); out$ci_lower_d <- q[1]; out$ci_upper_d <- q[2]
    } else {
      q <- score_interval(crn); out$ci_lower_b <- q[1]; out$ci_upper_b <- q[2]
    }
  }
  out$d_check_mean_z <- mean(check_z, na.rm = TRUE); out$d_check_n <- sum(is.finite(check_z))
  out$exact_median <- if (fy == "gaussian") median(exact, na.rm = TRUE) else NA_real_
  RNGkind(worker_kind[1L], worker_kind[2L], worker_kind[3L])
  assign(".Random.seed", worker_seed, envir = .GlobalEnv)
  out
}
RNGkind("L'Ecuyer-CMRG"); set.seed(base + task)
res <- parallel::mclapply(reps, one, mc.cores = cores)
out <- do.call(rbind, lapply(res, as.data.frame))
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_file, row.names = FALSE)
cat(sprintf("I6 task %d cell %s N=%d chunk=%d reps=%d boot=%d preflight=%s\n",
            task, sn, N, chunk, length(reps), boot_n, preflight))
