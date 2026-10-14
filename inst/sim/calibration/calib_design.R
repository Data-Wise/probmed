# Calibration design for the probmed single-mediator repair (tasks T15-T17).
#
# Implements the prespecified design in pmed's
# docs/SPEC-probmed-repair-2026-10-09.md ("Prespecified calibration design and
# decision rule"). Sourced by calib_fixtures.R, calib_worker.R,
# calib_benchmark.R and calib_collate.R. Truths are computed here from the
# generating parameters only; nothing in this file calls probmed.

CALIB_REPS <- 2000L
CALIB_N <- c(50L, 250L, 1000L)
# T22 (prespecified rerun after T21) sets CALIB_BASE_SEED=20261015 in the
# environment; the original run used the default.
CALIB_BASE_SEED <- as.integer(Sys.getenv("CALIB_BASE_SEED", "20261014"))

# --- MBCO matrix: 7 path pairs x 3 N = 21 cells ------------------------------
# X ~ Bernoulli(0.5), direct effect 0.3, Vm = 1, Vy = 1 (Vy = 2 for pair 7).
MBCO_PAIRS <- list(
  c(a = 0, b = 0), c(a = 0, b = 0.4), c(a = 0.5, b = 0), c(a = 0.05, b = 0.05),
  c(a = 0.5, b = 0.4), c(a = -0.5, b = 0.4), c(a = 2.5, b = 1.5)
)
mbco_vy <- function(k) if (k == 7L) 2 else 1

# --- Response-bootstrap matrix: 5 configurations x 3 N = 15 cells -------------
# X, C ~ Bernoulli(0.5) independent; Z ~ N(0, 1). Mediator predictor
# 0.5 X + 0.4 C (intercept 0 Gaussian, -0.25 Bernoulli logit). Outcome
# predictor intercept + 0.3 X + 0.4 M + 0.4 C (+ 0.3 Z in config 5).
# Bernoulli outcome intercepts are calibrated by quadrature (calib_fixtures.R)
# to marginal prevalence 0.05 (config 2) and 0.50 (configs 3, 5).
BOOT_CONFIGS <- list(
  list(id = 1L, m = "gaussian", y = "gaussian", z = FALSE, prevalence = NA),
  list(id = 2L, m = "gaussian", y = "binomial", z = FALSE, prevalence = 0.05),
  list(id = 3L, m = "gaussian", y = "binomial", z = FALSE, prevalence = 0.50),
  list(id = 4L, m = "binomial", y = "gaussian", z = FALSE, prevalence = NA),
  list(id = 5L, m = "binomial", y = "binomial", z = TRUE, prevalence = 0.50)
)
BOOT_N_BOOT <- 1000L
BOOT_N_SIM <- 10000L

calib_cells <- function() {
  mb <- expand.grid(k = seq_along(MBCO_PAIRS), N = CALIB_N)
  bt <- expand.grid(k = seq_along(BOOT_CONFIGS), N = CALIB_N)
  rbind(
    data.frame(cell_id = sprintf("mbco-%d-N%d", mb$k, mb$N), family = "mbco",
               k = mb$k, N = mb$N, stringsAsFactors = FALSE),
    data.frame(cell_id = sprintf("boot-%d-N%d", bt$k, bt$N), family = "boot",
               k = bt$k, N = bt$N, stringsAsFactors = FALSE)
  )
}

# --- RNG streams --------------------------------------------------------------
# L'Ecuyer-CMRG. Global replication index g = (cell_index - 1) * 2000 + rep;
# replication g uses stream g, i.e. nextRNGStream applied g times to the base
# stream set by set.seed(CALIB_BASE_SEED). Each replication starts from its own
# stream and pmed() is called with seed = NULL, so RNG consumption inside a
# replication never shifts another replication's stream.
calib_streams <- function(from, to) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(CALIB_BASE_SEED)
  s <- .Random.seed
  out <- vector("list", to - from + 1L)
  for (g in seq_len(to)) {
    s <- parallel::nextRNGStream(s)
    if (g >= from) out[[g - from + 1L]] <- s
  }
  out
}

# --- Data generation ----------------------------------------------------------
gen_mbco <- function(k, N) {
  ab <- MBCO_PAIRS[[k]]
  X <- stats::rbinom(N, 1, 0.5)
  M <- ab[["a"]] * X + stats::rnorm(N, 0, 1)
  Y <- 0.3 * X + ab[["b"]] * M + stats::rnorm(N, 0, sqrt(mbco_vy(k)))
  data.frame(X = X, M = M, Y = Y)
}

gen_boot <- function(k, N, fixtures) {
  cfg <- BOOT_CONFIGS[[k]]
  X <- stats::rbinom(N, 1, 0.5)
  C <- stats::rbinom(N, 1, 0.5)
  Z <- stats::rnorm(N)
  M <- if (cfg$m == "gaussian") {
    0.5 * X + 0.4 * C + stats::rnorm(N)
  } else {
    stats::rbinom(N, 1, stats::plogis(-0.25 + 0.5 * X + 0.4 * C))
  }
  iy <- boot_intercept_y(k, fixtures)
  eta <- iy + 0.3 * X + 0.4 * M + 0.4 * C + (if (cfg$z) 0.3 * Z else 0)
  Y <- if (cfg$y == "gaussian") eta + stats::rnorm(N) else stats::rbinom(N, 1, stats::plogis(eta))
  data.frame(X = X, C = C, Z = Z, M = M, Y = Y)
}

boot_intercept_y <- function(k, fixtures) {
  if (BOOT_CONFIGS[[k]]$y == "gaussian") return(0)
  fixtures$intercept_y[[paste0("boot-", k)]]
}

# --- Exact truths -------------------------------------------------------------
truth_mbco <- function(k) {
  ab <- MBCO_PAIRS[[k]]
  stats::pnorm(ab[["a"]] * ab[["b"]] / sqrt(2 * (ab[["b"]]^2 * 1 + mbco_vy(k))))
}

# Gaussian-mediator integral E[expit(lin + 0.4 M)], M ~ N(mu, 1).
.gm_expit <- function(lin, mu, rel.tol = 1e-12) {
  stats::integrate(function(m) stats::plogis(lin + 0.4 * m) * stats::dnorm(m, mu, 1),
                   -Inf, Inf, rel.tol = rel.tol)$value
}

# Truth standardized over this replication's observed C (and Z) rows.
truth_boot <- function(k, d, fixtures, rel.tol = 1e-12) {
  cfg <- BOOT_CONFIGS[[k]]
  if (cfg$m == "gaussian" && cfg$y == "gaussian") {
    # Additive Gaussian: covariates cancel; Y_t - Y_c ~ N(0.4*0.5, 2(0.16 + 1)).
    return(stats::pnorm(0.4 * 0.5 / sqrt(2 * (0.4^2 + 1))))
  }
  iy <- boot_intercept_y(k, fixtures)
  if (cfg$m == "binomial" && cfg$y == "gaussian") {
    # Y_t - Y_c = 0.4 (M_t - M_c) + N(0, 2); continuous, no ties.
    per_c <- vapply(c(0, 1), function(cc) {
      p1 <- stats::plogis(-0.25 + 0.5 + 0.4 * cc)
      p0 <- stats::plogis(-0.25 + 0.4 * cc)
      pd <- c(`1` = p1 * (1 - p0), `0` = p1 * p0 + (1 - p1) * (1 - p0),
              `-1` = (1 - p1) * p0)
      sum(pd * stats::pnorm(0.4 * c(1, 0, -1) / sqrt(2)))
    }, numeric(1))
    return(mean(per_c[d$C + 1]))
  }
  # Bernoulli outcome: half-tie P_med = 0.5 + 0.5 mean_r(p_t(r) - p_c(r)),
  # both arms evaluated with treatment 1 in the outcome and the same row r.
  arm_mean <- function(x_m, cc, z) {
    lin <- iy + 0.3 * 1 + 0.4 * cc + (if (cfg$z) 0.3 * z else 0)
    if (cfg$m == "gaussian") {
      .gm_expit(lin, 0.5 * x_m + 0.4 * cc, rel.tol)
    } else {
      pm <- stats::plogis(-0.25 + 0.5 * x_m + 0.4 * cc)
      pm * stats::plogis(lin + 0.4) + (1 - pm) * stats::plogis(lin)
    }
  }
  if (!cfg$z) {
    per_c <- vapply(c(0, 1), function(cc) arm_mean(1, cc, 0) - arm_mean(0, cc, 0), numeric(1))
    return(0.5 + 0.5 * mean(per_c[d$C + 1]))
  }
  diffs <- mapply(function(cc, z) arm_mean(1, cc, z) - arm_mean(0, cc, z), d$C, d$Z)
  0.5 + 0.5 * mean(diffs)
}

# Marginal outcome prevalence under the DGP for a candidate intercept.
boot_prevalence <- function(k, iy, rel.tol = 1e-12) {
  cfg <- BOOT_CONFIGS[[k]]
  cells <- expand.grid(x = 0:1, cc = 0:1)
  vals <- mapply(function(x, cc) {
    lin <- iy + 0.3 * x + 0.4 * cc
    if (cfg$m == "gaussian") {
      .gm_expit(lin, 0.5 * x + 0.4 * cc, rel.tol)
    } else {
      pm <- stats::plogis(-0.25 + 0.5 * x + 0.4 * cc)
      inner <- function(m) {
        if (!cfg$z) return(stats::plogis(lin + 0.4 * m))
        stats::integrate(function(z) stats::plogis(lin + 0.4 * m + 0.3 * z) * stats::dnorm(z),
                         -Inf, Inf, rel.tol = rel.tol)$value
      }
      pm * inner(1) + (1 - pm) * inner(0)
    }
  }, cells$x, cells$cc)
  mean(vals)
}

load_calib_fixtures <- function(path) {
  f <- read.dcf(path)
  ids <- grep("^intercept_y_boot_", colnames(f), value = TRUE)
  list(intercept_y = stats::setNames(
    as.list(as.numeric(f[1, ids])), gsub("_", "-", sub("^intercept_y_", "", ids))
  ))
}
