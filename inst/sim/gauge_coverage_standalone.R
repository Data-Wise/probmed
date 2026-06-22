#!/usr/bin/env Rscript
# STANDALONE coverage study for the gauge residual W and P_med (HPC-portable).
#
# Why standalone: the gauge estimator is pure base-`stats` (glm/lm); the S7
# class in probmed is only a container. To avoid installing the S7 + medfit +
# probmed chain on an HPC, this script INLINES `.corner_fit` (verbatim from
# R/corner.R) and the W / P_med + analytic-Wald / bootstrap-percentile interval
# math (verbatim from R/gauge-pmed.R, reps = 1). It must stay numerically
# identical to ward_residual() at the same seed -- verified by
# inst/sim/check_standalone_equiv.R before each run. Needs only stock R.
#
# Usage / output: identical CLI + columns to inst/sim/gauge_coverage.R.
#   Rscript gauge_coverage_standalone.R --smoke
#   SLURM array: one grid cell per $SLURM_ARRAY_TASK_ID -> $GAUGE_SIM_OUT/cell_<id>.rds

## ---- inlined estimator engine (verbatim from R/corner.R) -------------------
.corner_fit <- function(d, K, binY, covars) {
  cf <- paste(covars, collapse = " + ")
  f_pi <- stats::as.formula(paste("A ~", cf))
  f_q  <- stats::as.formula(paste("A ~ M +", cf))
  f_om <- stats::as.formula(paste("Y ~ A * M +", cf))
  n <- nrow(d); folds <- sample(rep(1:K, length.out = n))
  nm <- c("11", "10", "01", "00"); cor <- list(c(1, 1), c(1, 0), c(0, 1), c(0, 0))
  phi <- matrix(0, n, 4, dimnames = list(NULL, nm))
  gC <- numeric(n)
  for (k in 1:K) {
    tr <- d[folds != k, , drop = FALSE]
    te_i <- which(folds == k); te <- d[te_i, , drop = FALSE]
    pim <- stats::glm(f_pi, data = tr, family = stats::binomial())
    qm  <- stats::glm(f_q,  data = tr, family = stats::binomial())
    om  <- stats::glm(f_om, data = tr,
                      family = if (binY) stats::binomial() else stats::gaussian())
    p1 <- stats::predict(pim, newdata = te, type = "response"); gC[te_i] <- p1
    pa <- function(z) ifelse(z == 1, p1, 1 - p1)
    q1 <- stats::predict(qm, newdata = te, type = "response")
    qa <- function(z) ifelse(z == 1, q1, 1 - q1)
    mu <- function(z) stats::predict(om, newdata = transform(te, A = z), type = "response")
    for (j in 1:4) {
      a <- cor[[j]][1]; ap <- cor[[j]][2]
      muAM <- mu(a)
      ratio <- (qa(ap) / qa(a)) * (pa(a) / pa(ap))
      muAM_tr <- stats::predict(om, newdata = transform(tr, A = a), type = "response")
      sub <- tr$A == ap
      etam <- stats::lm(stats::reformulate(covars, "yy"),
                        data = cbind(data.frame(yy = muAM_tr[sub]),
                                     tr[sub, covars, drop = FALSE]))
      eta <- stats::predict(etam, newdata = te)
      phi[te_i, j] <- (te$A == a) / pa(a) * ratio * (te$Y - muAM) +
                      (te$A == ap) / pa(ap) * (muAM - eta) + eta
    }
  }
  list(phi = phi, g = gC)
}

## W, P_med and their intervals (verbatim logic from R/gauge-pmed.R, reps = 1).
gauge_core <- function(d, covars = "C", K = 5L, ci_level = 0.95,
                       se_method = c("analytic", "bootstrap"), B = 200L,
                       seed = 1L) {
  se_method <- match.arg(se_method)
  set.seed(seed)        # mirror ward_residual: re-seed before fold sampling
  binY <- all(d$Y %in% 0:1); n <- nrow(d)
  .gp_WP <- function(p) {
    t <- colMeans(p); oe <- t["11"] - t["00"]
    c(W = unname((oe - (t["10"] - t["00"]) - (t["01"] - t["00"])) / oe),
      P = unname((t["01"] - t["00"]) / oe))
  }
  phi <- .corner_fit(d, K, binY, covars)$phi
  th <- colMeans(phi)
  OE <- th["11"] - th["00"]; IDE <- th["10"] - th["00"]
  IIE <- th["01"] - th["00"]; R <- OE - IDE - IIE
  Pmed <- IIE / OE; W <- R / OE
  pOE <- phi[, "11"] - phi[, "00"]; pIDE <- phi[, "10"] - phi[, "00"]
  pIIE <- phi[, "01"] - phi[, "00"]; pR <- pOE - pIDE - pIIE
  se <- function(x) stats::sd(x) / sqrt(n); zc <- stats::qnorm(1 - (1 - ci_level) / 2)
  alpha <- 1 - ci_level
  seP <- se((pIIE - Pmed * pOE) / OE); seW <- se((pR - W * pOE) / OE)
  W_ci     <- c(W - zc * seW, W + zc * seW)
  p_med_ci <- c(Pmed - zc * seP, Pmed + zc * seP)
  if (se_method == "bootstrap") {
    bsamp <- vapply(seq_len(B), function(b) {
      db <- d[sample.int(n, n, replace = TRUE), , drop = FALSE]
      .gp_WP(.corner_fit(db, K, binY, covars)$phi)
    }, numeric(2))
    W_ci     <- stats::quantile(bsamp["W", ], c(alpha / 2, 1 - alpha / 2), names = FALSE)
    p_med_ci <- stats::quantile(bsamp["P", ], c(alpha / 2, 1 - alpha / 2), names = FALSE)
  }
  list(W = unname(W), P = unname(Pmed),
       W_ci = unname(W_ci), p_med_ci = unname(p_med_ci))
}

## ---- DGP + closed-form truth (same as gauge_coverage.R) --------------------
.pars <- list(g0 = -0.2, gC = 0.8, aM = 0.6, aC = 0.4,
              bA = 0.5, bM = 0.7, bC = 0.3)
true_WP <- function(tau, p = .pars) {
  OE  <- p$bA + (p$bM + tau) * p$aM
  c(W = (tau * p$aM) / OE, P = (p$bM * p$aM) / OE)
}
gen <- function(n, tau, seed, p = .pars) {
  set.seed(seed)
  C <- rnorm(n)
  A <- rbinom(n, 1, plogis(p$g0 + p$gC * C))
  M <- p$aM * A + p$aC * C + rnorm(n)
  Y <- p$bA * A + p$bM * M + tau * A * M + p$bC * C + rnorm(n)
  data.frame(A, M, Y, C)
}
one_rep <- function(n, tau, se_method, B, seed) {
  tr <- true_WP(tau)
  fit <- tryCatch(gauge_core(gen(n, tau, seed), se_method = se_method, B = B, seed = seed),
                  error = function(e) NULL)
  if (is.null(fit)) return(c(W = NA, P = NA))
  c(W = as.numeric(fit$W_ci[1] <= tr["W"] && tr["W"] <= fit$W_ci[2]),
    P = as.numeric(fit$p_med_ci[1] <= tr["P"] && tr["P"] <= fit$p_med_ci[2]))
}
grid <- function() {
  expand.grid(n = c(500L, 1000L, 2000L, 4000L),
              tau = c(0.0, 0.2, 0.8),
              se_method = c("analytic", "bootstrap"),
              stringsAsFactors = FALSE)
}
run_cell <- function(cell, nsim, B, base_seed = 1000L) {
  reps <- vapply(seq_len(nsim),
                 function(s) one_rep(cell$n, cell$tau, cell$se_method, B,
                                     seed = base_seed + s),
                 numeric(2))
  covW <- mean(reps["W", ], na.rm = TRUE); covP <- mean(reps["P", ], na.rm = TRUE)
  mcse <- function(p, m) sqrt(p * (1 - p) / m); tr <- true_WP(cell$tau)
  data.frame(n = cell$n, tau = cell$tau, se_method = cell$se_method,
             true_W = unname(tr["W"]), true_P = unname(tr["P"]),
             cov_W = covW, mcse_W = mcse(covW, nsim),
             cov_P = covP, mcse_P = mcse(covP, nsim), nsim = nsim, B = B)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  smoke <- "--smoke" %in% args
  getarg <- function(flag, default) {
    i <- match(flag, args); if (is.na(i)) default else as.integer(args[i + 1L])
  }
  nsim <- if (smoke) 30L else getarg("--nsim", 1000L)
  B    <- if (smoke) 100L else getarg("--B", 400L)
  g <- grid()
  task <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
  if (nzchar(task)) {
    id <- as.integer(task)
    outdir <- Sys.getenv("GAUGE_SIM_OUT", "gauge_sim_out")
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    res <- run_cell(g[id, ], nsim, B)
    saveRDS(res, file.path(outdir, sprintf("cell_%03d.rds", id)))
    message("wrote cell ", id); return(invisible(res))
  }
  if (smoke) g <- g[g$n == 500L & g$tau %in% c(0.0, 0.8), ]
  out <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) run_cell(g[i, ], nsim, B)))
  print(out, row.names = FALSE); invisible(out)
}
if (sys.nframe() == 0L) main()
