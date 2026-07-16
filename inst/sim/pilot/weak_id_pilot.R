## Pilot: how variable is weak_id_ratio across draws, and how often does the
## weak_id flag actually fire?  (Cited in ?GaugePmedResult.)
##
## PROVENANCE / STATUS: this is a 20-draw-per-point PILOT, not a measurement.
## At 20 draws a "50%" firing rate has a 95% interval of about [0.28, 0.72], and
## 0/20 alarms is only consistent with a false-alarm rate below ~15% -- NOT with
## zero. Committed so the numbers quoted in the roxygen are reproducible and so
## their weakness is auditable. The real study is
## ../hopper/run_weakid_validation.R (nsim = 2000/cell). Do not cite this file
## as validation.
##
## Run:  Rscript inst/sim/pilot/weak_id_pilot.R     (120 fits; ~25 min)
suppressMessages(devtools::load_all(quiet = TRUE))
expit <- function(x) 1 / (1 + exp(-x))

gen <- function(n, s, seed) {
  set.seed(seed)
  C <- rnorm(n); A <- rbinom(n, 1, expit(-0.2 + 0.8 * C))
  M <- s * 0.6 * A + 0.4 * C + rnorm(n)
  Y <- s * 0.5 * A + 0.7 * M + s * 0.4 * A * M + 0.3 * C + rnorm(n)
  data.frame(A, M, Y, C)
}

R <- 20L
rec <- list(); k <- 0L
cat(sprintf("%5s %8s %8s %9s %9s %8s\n",
            "s", "snr_mn", "snr_sd", "ratio_mn", "ratio_sd", "pct_flag"))
for (s in c(0.05, 0.10, 0.20, 0.30, 0.50, 1.00)) {
  snr <- rat <- flag <- numeric(R)
  for (i in seq_len(R)) {
    r <- suppressWarnings(
      ward_residual(gen(1500, s, 1000 + i), se_method = "bootstrap", B = 200L))
    snr[i] <- r@oe_snr; rat[i] <- r@weak_id_ratio; flag[i] <- isTRUE(r@weak_id)
    k <- k + 1L
    rec[[k]] <- data.frame(s = s, oe_snr = snr[i], ratio = rat[i], flag = flag[i])
  }
  cat(sprintf("%5.2f %8.2f %8.2f %9.2f %9.2f %7.0f%%\n",
              s, mean(snr), sd(snr), mean(rat), sd(rat), 100 * mean(flag)))
}
d <- do.call(rbind, rec)

## ---- PER-DRAW operating characteristics -------------------------------------
## The docs condition on oe_snr PER DRAW, not per s. Compute exactly that, with
## exact (Clopper-Pearson) intervals -- an earlier version of this script only
## printed per-s means, so the numbers quoted in ?GaugePmedResult were not
## actually the ones it computed. Pooling across s also means these rest on more
## draws than 20, so the intervals are tighter than a per-row reading implies.
ci <- function(k, n) if (n == 0) c(NA, NA) else stats::binom.test(k, n)$conf.int
weak   <- d$oe_snr <= 1.2
strong <- d$oe_snr >= 3
kw <- sum(d$flag[weak]);   nw <- sum(weak)
ks <- sum(d$flag[strong]); ns <- sum(strong)
cw <- ci(kw, nw); cs <- ci(ks, ns)
cat("\n---- per-draw OC (what the docs actually claim) ----\n")
cat(sprintf("sensitivity (flag | oe_snr <= 1.2): %d/%d = %.2f   exact 95%% CI [%.2f, %.2f]\n",
            kw, nw, kw / nw, cw[1], cw[2]))
cat(sprintf("false alarms (flag | oe_snr >= 3) : %d/%d = %.2f   exact 95%% CI [%.2f, %.2f]\n",
            ks, ns, ks / ns, cs[1], cs[2]))
cat(sprintf("ungrouped (1.2 < oe_snr < 3)      : %d draws\n", sum(!weak & !strong)))
cat("\nNOTE: ratio_sd ~= ratio_mn in the weak regime -- that variability is the\n",
    "mechanism behind the low sensitivity, and is better established than the\n",
    "sensitivity number itself. Cite the CIs above, not the point estimates.\n", sep = "")
