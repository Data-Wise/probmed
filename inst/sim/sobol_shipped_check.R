#!/usr/bin/env Rscript
# Post-weight-fix check of sobol_pmed() against exact truth: analytic se, the four
# interval constructions (Wald / A gated / B1 / B2), the boundary test, and (on a
# subset) the bootstrap se. One cell per call.
# Run from the package root. Usage: Rscript inst/sim/sobol_shipped_check.R <s> <tau> <n> [nrep=250] [B=100] [nb_rep=100]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
src <- readLines("inst/sim/phi_decomposition.R")
i1 <- grep("^gen <- function", src); i2 <- grep("^## ---- cells", src)
# eval(parse()) sources function DEFINITIONS (gen/truth/oracle) from a committed repo script.
eval(parse(text = src[i1:(i2 - 1)]))
args <- commandArgs(trailingOnly = TRUE)
s <- as.numeric(args[1]); tau <- as.numeric(args[2]); n <- as.integer(args[3])
nrep <- if (length(args) >= 4) as.integer(args[4]) else 250L
B <- if (length(args) >= 5) as.integer(args[5]) else 100L
nb_rep <- if (length(args) >= 6) as.integer(args[6]) else 100L
pd <- pm <- 0.5; cd <- pd * (1 - pd); cm <- pm * (1 - pm); cdm <- cd * cm
sob <- function(th) {                      # theta -> (Dm, Vm, VT, P), same algebra as .sobol_fit
  Dd <- (1 - pm) * (th["10"] - th["00"]) + pm * (th["11"] - th["01"])
  Dm <- (1 - pd) * (th["01"] - th["00"]) + pd * (th["11"] - th["10"])
  Rg <- th["11"] - th["10"] - th["01"] + th["00"]
  Vd <- cd * Dd^2; Vm <- cm * Dm^2; Vdm <- cdm * Rg^2; VT <- Vd + Vm + Vdm
  c(Dd = unname(Dd), Dm = unname(Dm), Rg = unname(Rg), Vm = unname(Vm), VT = unname(VT), P = unname(Vm / VT))
}
sob_if <- function(phi) {                  # analytic IF se for P and Dm from any influence matrix
  th <- colMeans(phi); o <- sob(th); cphi <- sweep(phi, 2, th)
  pDd <- (1 - pm) * (cphi[, "10"] - cphi[, "00"]) + pm * (cphi[, "11"] - cphi[, "01"])
  pDm <- (1 - pd) * (cphi[, "01"] - cphi[, "00"]) + pd * (cphi[, "11"] - cphi[, "10"])
  pR  <- cphi[, "11"] - cphi[, "10"] - cphi[, "01"] + cphi[, "00"]
  pVT <- 2 * cd * o["Dd"] * pDd + 2 * cm * o["Dm"] * pDm + 2 * cdm * o["Rg"] * pR
  pP  <- (2 * cm * o["Dm"] * pDm - o["P"] * pVT) / o["VT"]
  c(P = unname(o["P"]), se = sd(pP) / sqrt(nrow(phi)), Dm = unname(o["Dm"]), se_Dm = sd(pDm) / sqrt(nrow(phi)))
}
tr <- truth(s, tau, FALSE); truthS <- sob(tr$theta); Pt <- truthS["P"]
cover <- function(ci) unname(ci[1] <= Pt & Pt <= ci[2])   # unname: a named logical would mangle the vapply column names
t0 <- proc.time()[["elapsed"]]
rows <- t(vapply(seq_len(nrep), function(r) {
  seed <- 1000000L + r; set.seed(seed); d <- gen(n, s, tau, FALSE)
  f <- .sobol_fit(d, pd, pm, "C", 5L, seed = seed * 100L, warn_boundary = FALSE,
                  boundary_test = "split", procedure = "B", se_method = "analytic")
  so <- sob_if(phi_oracle(d, s, tau, FALSE))
  out <- c(P = f$P_med_sobol, se = f$se, Dm = f$Dm, se_Dm = f$se_Dm,
           cov_wald = cover(f$ci_wald), cov_A = cover(f$ci_A), cov_B1 = cover(f$ci_B1),
           cov_B2 = if (all(is.finite(f$ci_B2))) cover(f$ci_B2) else NA,
           w_wald = diff(f$ci_wald), w_B1 = diff(f$ci_B1), w_B2 = diff(f$ci_B2),
           boundary = as.numeric(f$boundary), reject = f$vmed_split_reject,
           P_or = unname(so["P"]), se_or = unname(so["se"]),
           cov_or = abs(so["P"] - Pt) <= qnorm(0.975) * so["se"],
           se_b = NA, cov_wald_b = NA, cov_B1_b = NA)
  if (r <= nb_rep && B > 0) {
    fb <- .sobol_fit(d, pd, pm, "C", 5L, seed = seed * 100L, warn_boundary = FALSE,
                     boundary_test = "split", procedure = "B", se_method = "bootstrap", B = B)
    out["se_b"] <- fb$se; out["cov_wald_b"] <- cover(fb$ci_wald); out["cov_B1_b"] <- cover(fb$ci_B1)
  }
  if (r %% 50 == 0) message(sprintf("  %d/%d  %.0fs", r, nrep, proc.time()[["elapsed"]] - t0))
  out
}, numeric(19)))
x <- as.data.frame(rows); emp <- sd(x$P); xb <- x[seq_len(min(nb_rep, nrep)), ]
cat(sprintf("sobol cell: s=%g tau=%g n=%d | nrep=%d B=%d (%d datasets) %.0fs | P_true=%.4f Dm_true=%.4f empSD(P)=%.4f sd(P_or)=%.4f mean(P)=%.4f\n",
            s, tau, n, nrep, B, nb_rep, proc.time()[["elapsed"]] - t0, Pt, truthS["Dm"], emp, sd(x$P_or), mean(x$P)))
out <- c(bias_P = mean(x$P) - Pt, se_over_empSD_med = median(x$se) / emp, se_cv = sd(x$se) / mean(x$se),
         oracle_se_over_empSD = median(x$se_or) / emp, oracle_cov = mean(x$cov_or),
         cov_wald = mean(x$cov_wald), cov_A = mean(x$cov_A), cov_B1 = mean(x$cov_B1), cov_B2 = mean(x$cov_B2, na.rm = TRUE),
         width_B1_over_wald_med = median(x$w_B1 / x$w_wald), width_B2_over_B1_med = median(x$w_B2 / x$w_B1, na.rm = TRUE),
         boundary_rate = mean(x$boundary), reject_rate = mean(x$reject),
         cov_wald_boundary = if (any(x$boundary == 1)) mean(x$cov_wald[x$boundary == 1]) else NA,
         cov_B1_boundary = if (any(x$boundary == 1)) mean(x$cov_B1[x$boundary == 1]) else NA,
         boot_se_over_empSD_med = median(xb$se_b, na.rm = TRUE) / emp, boot_se_cv = sd(xb$se_b, na.rm = TRUE) / mean(xb$se_b, na.rm = TRUE),
         cov_wald_boot = mean(xb$cov_wald_b, na.rm = TRUE), cov_B1_boot = mean(xb$cov_B1_b, na.rm = TRUE))
print(round(out, 3))
od <- Sys.getenv("GAUGE_SIM_OUT", "gauge_sim_out"); dir.create(od, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(rows = x, summary = out, cell = list(s = s, tau = tau, n = n, nrep = nrep, B = B, nb_rep = nb_rep, P_true = unname(Pt), Dm_true = unname(truthS["Dm"]))),
        file.path(od, sprintf("sobol_s%g_t%g_n%d.rds", s, tau, n)))
