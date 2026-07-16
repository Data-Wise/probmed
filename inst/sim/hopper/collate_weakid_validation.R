## Collate the #11 threshold-validation grid -> per-cell operating characteristics
## + the threshold sweep that actually validates (or corrects) the defaults.
parts <- list.files("~/weakid_val/parts", pattern = "^part_.*\\.rds$", full.names = TRUE)
if (!length(parts)) stop("no parts found in ~/weakid_val/parts")
df <- do.call(rbind, lapply(parts, readRDS))
mcse <- function(p, n) sqrt(p * (1 - p) / n)

## --- GUARD: the flag must actually vary. If weak_id is all-NA (old package, or
## se_method not live) or all-one-value, the sweep below is meaningless. -------
if (all(is.na(df$weak_id)))
  stop("GUARD: weak_id is all NA -- stale probmed (pre-#23) or bootstrap no-op. DO NOT TRUST.")
cat(sprintf("[GUARD] weak_id fires in %.1f%% of %d reps; weak_id_ratio range [%.2f, %.2f]\n",
            100 * mean(df$weak_id, na.rm = TRUE), nrow(df),
            min(df$weak_id_ratio, na.rm = TRUE), max(df$weak_id_ratio, na.rm = TRUE)))

## ---------------- per-cell summary ----------------
agg <- do.call(rbind, by(df, df[c("cell", "n", "s", "binY")], function(g) {
  ns <- nrow(g)
  data.frame(n = g$n[1], s = g$s[1], binY = g$binY[1], nsim = ns,
    trW = round(g$trW[1], 3),
    oe_snr_mn   = round(mean(g$oe_snr), 2),
    ratio_mn    = round(mean(g$weak_id_ratio, na.rm = TRUE), 2),
    ratio_sd    = round(sd(g$weak_id_ratio, na.rm = TRUE), 2),
    pct_weak_id = round(mean(g$weak_id, na.rm = TRUE), 3),
    pct_oe_irreg= round(mean(!g$oe_regular, na.rm = TRUE), 3),
    covW_an     = round(mean(g$covW_an), 3),
    mcse_covW_an= round(mcse(mean(g$covW_an), ns), 4),
    covW_pct    = round(mean(g$covW_pct), 3),
    row.names = NULL) }))
agg <- agg[order(agg$binY, agg$n, agg$s), ]

## ---------------- THE VALIDATION: coverage conditional on the flag ----------
## The flag claims "W's CI is least trustworthy here". The CI users get by
## default is the Wald one. So a threshold is validated iff, at that threshold,
## unflagged reps have ~nominal Wald coverage and flagged reps are clearly worse.
sweep <- do.call(rbind, lapply(seq(1.5, 6, by = 0.25), function(t) {
  fl <- df$weak_id_ratio >= t
  fl <- fl & !is.na(fl)
  if (sum(fl) < 50 || sum(!fl) < 50) return(NULL)   # too lopsided to read
  data.frame(threshold = t,
    pct_flagged   = round(mean(fl), 3),
    covW_flagged  = round(mean(df$covW_an[fl]), 3),
    covW_unflagged= round(mean(df$covW_an[!fl]), 3),
    separation    = round(mean(df$covW_an[!fl]) - mean(df$covW_an[fl]), 3),
    row.names = NULL) }))

## Same for the regularity gate, on its own scale.
sweep_snr <- do.call(rbind, lapply(seq(1, 4, by = 0.25), function(t) {
  fl <- df$oe_snr < t
  if (sum(fl) < 50 || sum(!fl) < 50) return(NULL)
  data.frame(threshold = t,
    pct_flagged   = round(mean(fl), 3),
    covW_flagged  = round(mean(df$covW_an[fl]), 3),
    covW_unflagged= round(mean(df$covW_an[!fl]), 3),
    separation    = round(mean(df$covW_an[!fl]) - mean(df$covW_an[fl]), 3),
    row.names = NULL) }))

## ---------------- sensitivity / specificity at the SHIPPED defaults ---------
## "weak regime" = the cells the flag is meant to catch; "strong" = it must not.
weak   <- df$oe_snr <= 1.2
strong <- df$oe_snr >= 3
oc <- data.frame(
  sensitivity_at_3 = round(mean(df$weak_id[weak],   na.rm = TRUE), 3),  # want high
  fpr_at_3         = round(mean(df$weak_id[strong], na.rm = TRUE), 3),  # want ~0
  n_weak = sum(weak), n_strong = sum(strong))

dir.create("~/weakid_val", showWarnings = FALSE)
write.csv(agg,       "~/weakid_val/weakid_validation_cells.csv", row.names = FALSE)
write.csv(sweep,     "~/weakid_val/weakid_threshold_sweep.csv",  row.names = FALSE)
write.csv(sweep_snr, "~/weakid_val/oe_snr_threshold_sweep.csv",  row.names = FALSE)
saveRDS(df,          "~/weakid_val/weakid_validation_raw.rds")

cat("\n==== per-cell ====\n");                     print(agg)
cat("\n==== weak_id_ratio_threshold sweep ====\n");print(sweep)
cat("\n==== oe_snr_threshold sweep ====\n");       print(sweep_snr)
cat("\n==== operating characteristics @ shipped defaults (3 / 2) ====\n"); print(oc)
cat(sprintf("\ncollated %d reps across %d cells\n", nrow(df), nrow(agg)))
