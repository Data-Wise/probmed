#!/usr/bin/env Rscript
# Post-weight-fix check of pmedW_dr(): bias / empSD of the point estimate pre- vs
# post-fix against a simulated exact truth, and bootstrap coverage via
# wasserstein_pmed(method = "dr") on a subset. Continuous Y only.
# Run from the package root. Usage: Rscript inst/sim/pmedw_shipped_check.R <s> <tau> <n> [nrep=200] [n_boot=100] [nb_rep=60]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
src <- readLines("inst/sim/phi_decomposition.R")
i1 <- grep("^gen <- function", src); i2 <- grep("^## ---- cells", src)
# eval(parse()) sources function DEFINITIONS (gen) from a committed repo script.
eval(parse(text = src[i1:(i2 - 1)]))
## pre-fix pmedW_dr, sourced from the parent of the fix commit into its own env
old <- new.env()
old_src <- system2("git", c("-C", ".", "show", "a73efc8^:R/wasserstein-pmed.R"), stdout = TRUE)
eval(parse(text = old_src), envir = old)   # function DEFINITIONS from this repo's own history
args <- commandArgs(trailingOnly = TRUE)
s <- as.numeric(args[1]); tau <- as.numeric(args[2]); n <- as.integer(args[3])
nrep <- if (length(args) >= 4) as.integer(args[4]) else 200L
n_boot <- if (length(args) >= 5) as.integer(args[5]) else 100L
nb_rep <- if (length(args) >= 6) as.integer(args[6]) else 60L
## exact-ish truth: simulate the corner laws Y(a, M(a')) from the DGP equations
set.seed(424242); N <- 2e6
C <- rnorm(N); eM <- rnorm(N); eY <- rnorm(N)
corner <- function(a, ap) { M <- 0.6 * s * ap + 0.4 * C + eM; 0.5 * s * a + 0.7 * M + tau * s * a * M + 0.3 * C + eY }
nu11 <- corner(1, 1); nu10 <- corner(1, 0); nu00 <- corner(0, 0)
NIEt <- .w2_1d(nu11, nu10); NDEt <- .w2_1d(nu10, nu00); Pt <- NIEt / (NIEt + NDEt)
rm(C, eM, eY, nu11, nu10, nu00)
t0 <- proc.time()[["elapsed"]]
rows <- t(vapply(seq_len(nrep), function(r) {
  seed <- 1000000L + r; set.seed(seed); d <- gen(n, s, tau, FALSE)
  new <- pmedW_dr(d, covars = "C", K = 5L, G = 300L, seed = seed * 100L)
  pre <- old$pmedW_dr(d, covars = "C", K = 5L, G = 300L, seed = seed * 100L)
  out <- c(P_new = new$pmedW, NIE_new = new$NIE_W, NDE_new = new$NDE_W, P_pre = pre$pmedW,
           se_b = NA, lo_b = NA, hi_b = NA)
  if (r <= nb_rep && n_boot > 0) {
    w <- suppressWarnings(suppressMessages(wasserstein_pmed(d, covars = "C", method = "dr", n_boot = n_boot,
                                                            K = 5L, G = 300L, seed = seed * 100L)))
    out["se_b"] <- w@pmedW_se; out["lo_b"] <- w@pmedW_ci[1]; out["hi_b"] <- w@pmedW_ci[2]
  }
  if (r %% 25 == 0) message(sprintf("  %d/%d  %.0fs", r, nrep, proc.time()[["elapsed"]] - t0))
  out
}, numeric(7)))
x <- as.data.frame(rows); xb <- x[seq_len(min(nb_rep, nrep)), ]
cat(sprintf("pmedW cell: s=%g tau=%g n=%d | nrep=%d n_boot=%d (%d datasets) %.0fs | P_true=%.4f NIE=%.4f NDE=%.4f\n",
            s, tau, n, nrep, n_boot, nb_rep, proc.time()[["elapsed"]] - t0, Pt, NIEt, NDEt))
out <- c(bias_post = mean(x$P_new) - Pt, bias_pre = mean(x$P_pre) - Pt,
         empSD_post = sd(x$P_new), empSD_pre = sd(x$P_pre), rmse_post = sqrt(mean((x$P_new - Pt)^2)), rmse_pre = sqrt(mean((x$P_pre - Pt)^2)),
         boot_se_over_empSD_med = median(xb$se_b, na.rm = TRUE) / sd(x$P_new), boot_se_cv = sd(xb$se_b, na.rm = TRUE) / mean(xb$se_b, na.rm = TRUE),
         cov_boot = mean(xb$lo_b <= Pt & Pt <= xb$hi_b, na.rm = TRUE), n_boot_datasets = sum(!is.na(xb$se_b)))
print(round(out, 4))
od <- Sys.getenv("GAUGE_SIM_OUT", "gauge_sim_out"); dir.create(od, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(rows = x, summary = out, cell = list(s = s, tau = tau, n = n, nrep = nrep, n_boot = n_boot, nb_rep = nb_rep, P_true = Pt, NIE = NIEt, NDE = NDEt)),
        file.path(od, sprintf("pmedw_s%g_t%g_n%d.rds", s, tau, n)))
