# Collate sobol_shipped_check.R outputs. Usage: Rscript inst/sim/sobol_collate.R <dir with sobol_*.rds> [out.csv]
S <- commandArgs(trailingOnly = TRUE)[1]
fs <- list.files(S, "^sobol_.*[.]rds$", full.names = TRUE)
summ <- function(f) {
  o <- readRDS(f); x <- o$rows; names(x) <- sub("[.]P$", "", names(x)); cl <- o$cell
  Pt <- cl$P_true; emp <- sd(x$P); xb <- x[seq_len(min(cl$nb_rep, cl$nrep)), ]
  c(s = cl$s, tau = cl$tau, n = cl$n, P_true = round(Pt, 4), Dm_true = round(cl$Dm_true, 3),
    bias = round(mean(x$P) - Pt, 4), empSD = round(emp, 4), oracle_SD = round(sd(x$P_or), 4),
    se_ratio = round(median(x$se) / emp, 3), se_cv = round(sd(x$se) / mean(x$se), 3),
    or_se_ratio = round(median(x$se_or) / emp, 3), cov_oracle = mean(x$cov_or),
    cov_wald = mean(x$cov_wald), cov_A = mean(x$cov_A), cov_B1 = mean(x$cov_B1), cov_B2 = round(mean(x$cov_B2, na.rm = TRUE), 3),
    wB1_wald = round(median(x$w_B1 / x$w_wald), 3), wB2_B1 = round(median(x$w_B2 / x$w_B1, na.rm = TRUE), 3),
    boundary = mean(x$boundary), reject = mean(x$reject),
    cov_wald_bnd = if (any(x$boundary == 1)) round(mean(x$cov_wald[x$boundary == 1]), 3) else NA,
    cov_B1_bnd = if (any(x$boundary == 1)) round(mean(x$cov_B1[x$boundary == 1]), 3) else NA,
    cov_wald_nobnd = if (any(x$boundary == 0)) round(mean(x$cov_wald[x$boundary == 0]), 3) else NA,
    boot_se_ratio = round(median(xb$se_b, na.rm = TRUE) / emp, 3), boot_se_cv = round(sd(xb$se_b, na.rm = TRUE) / mean(xb$se_b, na.rm = TRUE), 3),
    cov_wald_boot = round(mean(xb$cov_wald_b, na.rm = TRUE), 3), cov_B1_boot = round(mean(xb$cov_B1_b, na.rm = TRUE), 3),
    nrep = cl$nrep, nb = sum(!is.na(xb$se_b)))
}
out <- as.data.frame(do.call(rbind, lapply(fs, summ))); out <- out[order(-out$s, out$tau, out$n), ]
print(t(out), quote = FALSE)
if (length(commandArgs(trailingOnly = TRUE)) >= 2) write.csv(out, commandArgs(trailingOnly = TRUE)[2], row.names = FALSE)
