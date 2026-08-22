#!/usr/bin/env Rscript
# SE candidates, round 2 (cell 1: n=800, cont, tau=0; reps=10 point estimator).
# Round 1 showed every partition-aggregation candidate has CV 0.54-0.66 vs the
# oracle's 0.18: the instability comes from dataset-level nuisance-fit quality,
# not from how partitions are combined. Round 2 targets that directly:
#   * IF sd from FULL-SAMPLE nuisance fits (less noisy than fold fits), used for
#     the variance only -- the point estimate stays the cross-fit reps=10 one;
#   * winsorized IF contributions (1%/99%), on both the full-sample and the
#     partition-averaged phi, in case a few extreme contributions drive the CV.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
src <- readLines("inst/sim/phi_decomposition.R")
i1 <- grep("^gen <- function", src); i2 <- grep("^## ---- cells", src)
# eval(parse()) sources function DEFINITIONS from a committed repo script, not
# external input -- the standard way to reuse a CLI script's functions.
eval(parse(text = src[i1:(i2 - 1)]))

args <- commandArgs(trailingOnly = TRUE)
n <- 800L; K <- 5L; R <- 10L
nrep <- if (length(args)) as.integer(args[1]) else 250L
s <- 1; tau <- 0; binY <- FALSE; covars <- "C"
tr <- truth(s, tau, binY)

Wof <- function(phi) { t <- colMeans(phi); OE <- t["11"] - t["00"]
  unname((OE - (t["10"] - t["00"]) - (t["01"] - t["00"])) / OE) }
ifW <- function(phi, W = NULL) { t <- colMeans(phi); OE <- t["11"] - t["00"]
  if (is.null(W)) W <- (OE - (t["10"] - t["00"]) - (t["01"] - t["00"])) / OE
  pOE <- phi[, "11"] - phi[, "00"]; pR <- pOE - (phi[, "10"] - phi[, "00"]) - (phi[, "01"] - phi[, "00"])
  (pR - W * pOE) / OE }
wins <- function(v, p = 0.01) { q <- quantile(v, c(p, 1 - p)); pmin(pmax(v, q[1]), q[2]) }
folds_of <- function(seed) { set.seed(seed); sample(rep(1:K, length.out = n)) }
clust_se <- function(psi, f) { psi <- psi - mean(psi); S <- tapply(psi, f, sum); sqrt(sum(S^2)) / n }

## full-sample nuisance fits: .corner_fit's inner loop with train = test = d
full_phi <- function(d) {
  cf <- paste(covars, collapse = " + ")
  pim <- glm(as.formula(paste("A ~", cf)), data = d, family = binomial())
  qm  <- glm(as.formula(paste("A ~ M +", cf)), data = d, family = binomial())
  om  <- glm(as.formula(paste("Y ~ A * M +", cf)), data = d, family = if (binY) binomial() else gaussian())
  p1 <- fitted(pim); q1 <- fitted(qm)
  pa <- function(z) if (z == 1) p1 else 1 - p1; qa <- function(z) if (z == 1) q1 else 1 - q1
  mu <- function(z) predict(om, newdata = transform(d, A = z), type = "response")
  nm <- c("11", "10", "01", "00"); cor <- list(c(1, 1), c(1, 0), c(0, 1), c(0, 0))
  phi <- matrix(0, n, 4, dimnames = list(NULL, nm))
  for (j in 1:4) { a <- cor[[j]][1]; ap <- cor[[j]][2]
    muAM <- mu(a); ratio <- (qa(ap) / qa(a)) * (pa(a) / pa(ap))
    sub <- d$A == ap
    etam <- lm(reformulate(covars, "yy"), data = cbind(data.frame(yy = muAM[sub]), d[sub, covars, drop = FALSE]))
    eta <- predict(etam, newdata = d)
    phi[, j] <- (d$A == a) / pa(a) * ratio * (d$Y - muAM) + (d$A == ap) / pa(ap) * (muAM - eta) + eta }
  phi
}

rows <- vector("list", nrep); t0 <- proc.time()[["elapsed"]]
for (r in seq_len(nrep)) {
  seed <- 1000000L + r; set.seed(seed); d <- gen(n, s, tau, binY)
  phis <- vector("list", R); fl <- vector("list", R)
  for (p in seq_len(R)) { ps <- seed * 100L + p; fl[[p]] <- folds_of(ps); set.seed(ps)
    phis[[p]] <- .corner_fit(d, K, binY, covars)$phi }
  pbar <- Reduce(`+`, phis) / R; W_avg <- Wof(pbar)
  pf <- full_phi(d)
  ifp <- lapply(seq_len(R), function(p) ifW(phis[[p]]))
  seCluM  <- median(sapply(seq_len(R), function(p) clust_se(ifp[[p]], fl[[p]])))   # round-1 best
  seFull  <- sd(ifW(pf, W_avg)) / sqrt(n)                 # full-sample nuisances, cross-fit W
  seFullW <- sd(wins(ifW(pf, W_avg))) / sqrt(n)           # + winsorized 1/99
  seBarW  <- sd(wins(ifW(pbar))) / sqrt(n)                # averaged phi, winsorized
  seCluMW <- median(sapply(seq_len(R), function(p) clust_se(wins(ifp[[p]]), fl[[p]])))
  so <- gauge_stats(phi_oracle(d, s, tau, binY)); ior <- ifW(phi_oracle(d, s, tau, binY))
  seOrW <- sd(wins(ior)) / sqrt(n)
  rows[[r]] <- c(W = W_avg, W_or = unname(so["W"]), se_or = unname(so["seW"]), se_orW = seOrW,
                 seCluM = seCluM, seFull = seFull, seFullW = seFullW, seBarW = seBarW, seCluMW = seCluMW)
}
x <- as.data.frame(do.call(rbind, rows)); emp <- sd(x$W); z <- qnorm(0.975)
score <- function(se, W = x$W) c(med_ratio = median(se) / emp, mean_ratio = mean(se) / emp,
  cv = sd(se) / mean(se), cov = mean(abs(W - tr$W) <= z * se),
  q05 = quantile(se, .05) / emp, q95 = quantile(se, .95) / emp)
res <- rbind(fold_clustered_median = score(x$seCluM),
             fullsample_IF = score(x$seFull), fullsample_IF_wins = score(x$seFullW),
             avgphi_IF_wins = score(x$seBarW), fold_clustered_median_wins = score(x$seCluMW),
             oracle = score(x$se_or, x$W_or), oracle_wins = score(x$se_orW, x$W_or))
cat(sprintf("round 2, cell 1, reps=%d, nrep=%d, %.0fs | empSD(W_avg)=%.4f sd(W_or)=%.4f\n",
            R, nrep, proc.time()[["elapsed"]] - t0, emp, sd(x$W_or)))
print(round(res, 3))
saveRDS(list(rows = x, res = res), file.path(Sys.getenv("GAUGE_SIM_OUT", "gauge_sim_out"), "se_candidates2_cell1.rds"))
