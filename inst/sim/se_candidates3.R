#!/usr/bin/env Rscript
# SE candidates, round 3: does the round-2 winner (IF variance from FULL-SAMPLE
# nuisance fits, cross-fit reps=10 point estimate) hold outside cell 1?
# Usage: Rscript se_candidates3.R <binY 0/1> <s> <tau> [nrep]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
src <- readLines("inst/sim/phi_decomposition.R")
i1 <- grep("^gen <- function", src); i2 <- grep("^## ---- cells", src)
# eval(parse()) sources function DEFINITIONS from a committed repo script, not
# external input -- the standard way to reuse a CLI script's functions.
eval(parse(text = src[i1:(i2 - 1)]))

args <- commandArgs(trailingOnly = TRUE)
binY <- as.integer(args[1]) == 1L; s <- as.numeric(args[2]); tau <- as.numeric(args[3])
nrep <- if (length(args) >= 4) as.integer(args[4]) else 250L
n <- 800L; K <- 5L; R <- 10L; covars <- "C"
tr <- truth(s, tau, binY)

Wof <- function(phi) { t <- colMeans(phi); OE <- t["11"] - t["00"]
  unname((OE - (t["10"] - t["00"]) - (t["01"] - t["00"])) / OE) }
ifW <- function(phi, W = NULL) { t <- colMeans(phi); OE <- t["11"] - t["00"]
  if (is.null(W)) W <- (OE - (t["10"] - t["00"]) - (t["01"] - t["00"])) / OE
  pOE <- phi[, "11"] - phi[, "00"]; pR <- pOE - (phi[, "10"] - phi[, "00"]) - (phi[, "01"] - phi[, "00"])
  (pR - W * pOE) / OE }
folds_of <- function(seed) { set.seed(seed); sample(rep(1:K, length.out = n)) }
clust_se <- function(psi, f) { psi <- psi - mean(psi); S <- tapply(psi, f, sum); sqrt(sum(S^2)) / n }

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
  pbar <- Reduce(`+`, phis) / R; W_avg <- Wof(pbar); Wr <- sapply(phis, Wof)
  seIF <- sd(ifW(pbar)) / sqrt(n); seMC <- sqrt(var(Wr) / R); seShip <- sqrt(seIF^2 + seMC^2)
  ifp <- lapply(seq_len(R), function(p) ifW(phis[[p]]))
  seCluM <- median(sapply(seq_len(R), function(p) clust_se(ifp[[p]], fl[[p]])))
  seFull <- sd(ifW(full_phi(d), W_avg)) / sqrt(n)
  so <- gauge_stats(phi_oracle(d, s, tau, binY))
  rows[[r]] <- c(W = W_avg, W_or = unname(so["W"]), se_or = unname(so["seW"]),
                 seShip = seShip, seCluM = seCluM, seFull = seFull)
}
x <- as.data.frame(do.call(rbind, rows)); emp <- sd(x$W); z <- qnorm(0.975)
score <- function(se, W = x$W) c(med_ratio = median(se) / emp, mean_ratio = mean(se) / emp,
  cv = sd(se) / mean(se), cov = mean(abs(W - tr$W) <= z * se),
  q05 = quantile(se, .05) / emp, q95 = quantile(se, .95) / emp)
res <- rbind(shipped_total = score(x$seShip), fold_clustered_median = score(x$seCluM),
             fullsample_IF = score(x$seFull), oracle = score(x$se_or, x$W_or))
cat(sprintf("round 3: n=%d binY=%s s=%g tau=%g, reps=%d, nrep=%d, %.0fs | W_true=%.4f empSD(W_avg)=%.4f sd(W_or)=%.4f\n",
            n, binY, s, tau, R, nrep, proc.time()[["elapsed"]] - t0, tr$W, emp, sd(x$W_or)))
print(round(res, 3))
saveRDS(list(rows = x, res = res), sprintf(
  file.path(Sys.getenv("GAUGE_SIM_OUT", "gauge_sim_out"), "se_candidates3_b%d_s%g_t%g.rds"), binY, s, tau))
