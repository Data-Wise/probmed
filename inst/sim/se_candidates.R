#!/usr/bin/env Rscript
# SE candidates for the reps > 1 gauge point estimator (cell 1: n=800, cont, tau=0).
# Target: the oracle's behavior -- median ratio ~1, CV ~0.17, Wald coverage ~0.95.
# Runs against the dev tree via pkgload; reuses the decomposition script's functions.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
src <- readLines("inst/sim/phi_decomposition.R")
i1 <- grep("^gen <- function", src); i2 <- grep("^## ---- cells", src)
# eval(parse()) here sources function DEFINITIONS (gen, GH, oracle nuisances,
# truth, gauge_stats) from a committed script in this repo, not external input --
# the standard way to reuse a CLI script's functions without exporting them.
eval(parse(text = src[i1:(i2 - 1)]))

args <- commandArgs(trailingOnly = TRUE)
n <- 800L; K <- 5L; R <- 10L
nrep <- if (length(args)) as.integer(args[1]) else 250L
s <- 1; tau <- 0; binY <- FALSE
tr <- truth(s, tau, binY)

Wof <- function(phi) { t <- colMeans(phi); OE <- t["11"] - t["00"]
  unname((OE - (t["10"] - t["00"]) - (t["01"] - t["00"])) / OE) }
ifW <- function(phi, W = NULL) { t <- colMeans(phi); OE <- t["11"] - t["00"]
  if (is.null(W)) W <- (OE - (t["10"] - t["00"]) - (t["01"] - t["00"])) / OE
  pOE <- phi[, "11"] - phi[, "00"]; pR <- pOE - (phi[, "10"] - phi[, "00"]) - (phi[, "01"] - phi[, "00"])
  (pR - W * pOE) / OE }
folds_of <- function(seed) { set.seed(seed); sample(rep(1:K, length.out = n)) }   # .corner_fit's first draw
clust_se <- function(psi, f) {                   # fold-clustered sandwich for a mean of psi
  psi <- psi - mean(psi); S <- tapply(psi, f, sum); sqrt(sum(S^2)) / n }

rows <- vector("list", nrep); t0 <- proc.time()[["elapsed"]]
for (r in seq_len(nrep)) {
  seed <- 1000000L + r; set.seed(seed); d <- gen(n, s, tau, binY)
  phis <- vector("list", R); fl <- vector("list", R)
  for (p in seq_len(R)) { ps <- seed * 100L + p; fl[[p]] <- folds_of(ps); set.seed(ps)
    phis[[p]] <- .corner_fit(d, K, binY, "C")$phi }
  pbar <- Reduce(`+`, phis) / R; W_avg <- Wof(pbar); Wr <- sapply(phis, Wof)
  ## candidates
  seIF   <- sd(ifW(pbar)) / sqrt(n)                                   # shipped, IF part
  seMC   <- sqrt(var(Wr) / R)                                         # shipped, partition add-on
  seShip <- sqrt(seIF^2 + seMC^2)                                     # shipped total
  ifp    <- lapply(seq_len(R), function(p) ifW(phis[[p]]))            # per-partition IF, own W
  se_p   <- sapply(ifp, sd) / sqrt(n)
  seMean <- mean(se_p); seMed <- median(se_p)                         # mean / median over partitions
  seOwn  <- sd(Reduce(`+`, ifp) / R) / sqrt(n)                        # per-obs average of own-W IFs
  seClu  <- mean(sapply(seq_len(R), function(p) clust_se(ifp[[p]], fl[[p]])))   # fold-clustered, averaged
  seCluM <- median(sapply(seq_len(R), function(p) clust_se(ifp[[p]], fl[[p]])))
  so <- gauge_stats(phi_oracle(d, s, tau, binY))
  rows[[r]] <- c(W = W_avg, W_or = unname(so["W"]), se_or = unname(so["seW"]),
                 seIF = seIF, seMC = seMC, seShip = seShip, seMean = seMean, seMed = seMed,
                 seOwn = seOwn, seClu = seClu, seCluM = seCluM)
}
x <- as.data.frame(do.call(rbind, rows)); emp <- sd(x$W); z <- qnorm(0.975)
score <- function(se, W = x$W) c(med_ratio = median(se) / emp, mean_ratio = mean(se) / emp,
  cv = sd(se) / mean(se), cov = mean(abs(W - tr$W) <= z * se),
  q05 = quantile(se, .05) / emp, q95 = quantile(se, .95) / emp)
res <- rbind(
  shipped_IF_part = score(x$seIF), shipped_MC_part = score(x$seMC), shipped_total = score(x$seShip),
  perpart_mean = score(x$seMean), perpart_median = score(x$seMed), perpart_ownW_avg = score(x$seOwn),
  fold_clustered_mean = score(x$seClu), fold_clustered_median = score(x$seCluM),
  oracle = score(x$se_or, x$W_or))
cat(sprintf("cell 1 (n=%d cont tau=0), reps=%d, nrep=%d, %.0fs | empSD(W_avg)=%.4f  sd(W_or)=%.4f\n",
            n, R, nrep, proc.time()[["elapsed"]] - t0, emp, sd(x$W_or)))
print(round(res, 3))
saveRDS(list(rows = x, res = res), file.path(Sys.getenv("GAUGE_SIM_OUT", "gauge_sim_out"), "se_candidates_cell1.rds"))
