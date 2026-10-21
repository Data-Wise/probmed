## Collate gauge bootstrap grid -> per-cell coverage with Monte Carlo SE (gate B1).
bdir <- Sys.getenv("GAUGE_BOOT_DIR", "~/gauge_boot")   # post-fix rerun: ~/gauge_boot_postfix
parts <- list.files(file.path(bdir, "parts"), pattern="^part_.*\\.rds$", full.names=TRUE)
df <- do.call(rbind, lapply(parts, readRDS))
mcse <- function(p,n) sqrt(p*(1-p)/n)
agg <- do.call(rbind, by(df, df[c("cell","n","tint","binY")], function(g){
  ns<-nrow(g)
  data.frame(n=g$n[1], tint=g$tint[1], binY=g$binY[1], nsim=ns,
    trW=round(g$trW[1],3), biasW=round(mean(g$W)-g$trW[1],3), empSD_W=round(sd(g$W),3),
    covW_an=round(mean(g$covW_an),3),  mcse_covW_an =round(mcse(mean(g$covW_an),ns),4),
    covW_pct=round(mean(g$covW_pct),3),mcse_covW_pct=round(mcse(mean(g$covW_pct),ns),4),
    covP_an=round(mean(g$covP_an),3),  covP_pct=round(mean(g$covP_pct),3),
    seW_an_ratio=round(mean(g$seW_an)/sd(g$W),2), seW_bt_ratio=round(mean(g$seW_bt)/sd(g$W),2),
    divW_med=round(median(g$divW_WaldVsPct),3), row.names=NULL) }))
## --- STALE-PACKAGE GUARD (WATCH-OUTS #1): the bootstrap arm MUST differ from analytic ---
## If se_method="bootstrap" silently no-ops (wrong R_LIBS / stale installed probmed), seW_bt == seW_an
## and the percentile coverage is really the OLD analytic CI. Flag unmistakably; never silently ship.
bt_eq_an    <- mean(abs(df$seW_bt - df$seW_an) < 1e-9, na.rm = TRUE)
bt_an_ratio <- median(df$seW_bt / df$seW_an, na.rm = TRUE)
cat(sprintf("[GUARD] bootstrap==analytic in %.1f%% of reps; median seW_bt/seW_an = %.2f\n",
            100 * bt_eq_an, bt_an_ratio))
## Post-fix (2026-08-22) the bootstrap se legitimately ~= the analytic se (ratio ~1.0-1.1;
## the pre-fix 1.82 was the weight bug), so only EXACT equality indicates a no-op arm.
if (bt_eq_an > 0.5 || !is.finite(bt_an_ratio)) {
  msg <- sprintf("STALE-PACKAGE SUSPECTED: bootstrap se ~= analytic (eq=%.1f%%, ratio=%.2f) -- se_method='bootstrap' may be a NO-OP (wrong R_LIBS). DO NOT TRUST coverage.",
                 100 * bt_eq_an, bt_an_ratio)
  cat("[GUARD] *** ", msg, " ***\n", sep = "")
  writeLines(msg, file.path(bdir, "STALE_WARNING.txt"))
} else {
  cat("[GUARD] OK -- bootstrap se genuinely differs from analytic (live se_method='bootstrap').\n")
}
write.csv(agg, file.path(bdir, "gauge_boot_grid.csv"), row.names=FALSE)
saveRDS(df,  file.path(bdir, "gauge_boot_raw.rds"))   # for SE-vs-estimate scatter (gate B2)
cat("collated", nrow(df), "reps across", nrow(agg), "cells -> gauge_boot_grid.csv\n"); print(agg)
