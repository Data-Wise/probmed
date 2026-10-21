## Truth functions for the RG-flow estimator tests.
##
## These are a self-contained reimplementation of `pop_paths()` / `C0_analytic()`
## from the manuscript repo
## (research/mediation-noncollapsibility, 01-rg-pmed/sims/dgp.R, verified
## kill-test 2026-06-23). probmed cannot depend on that repo, so the closed
## forms are vendored here; the hardcoded constants in test-rg-flow.R were
## generated from the ORIGINAL there and are cross-checked against this
## reimplementation, so a transcription slip in either copy fails the suite.
##
## Two-level linear-Gaussian DGP, parameters `p`:
##   tb, tw   between/within variance of A
##   ab, aw   a-path (A -> M) between/within
##   bb, bw   b-path (M -> Y) between/within
##   cb, cw   direct path between/within
##   sMb,sMw  mediator residual variance between/within
##   sYb,sYw  outcome residual variance (drops out of g)

## Population least-squares paths on cluster-mean data of block size m.
rg_pop_paths <- function(p, m) {
  w <- 1 / m
  VA <- p$tb + w * p$tw
  CAM <- p$ab * p$tb + w * p$aw * p$tw
  tbf <- p$ab * p$bb + p$cb
  twf <- p$aw * p$bw + p$cw
  CAY <- tbf * p$tb + w * twf * p$tw
  VM <- (p$ab^2 * p$tb + p$sMb) + w * (p$aw^2 * p$tw + p$sMw)
  CMY <- (p$bb * (p$ab^2 * p$tb + p$sMb) + p$cb * p$ab * p$tb) +
    w * (p$bw * (p$aw^2 * p$tw + p$sMw) + p$cw * p$aw * p$tw)
  a <- CAM / VA
  den <- VM * VA - CAM^2
  b <- (CMY * VA - CAY * CAM) / den
  TE <- CAY / VA
  list(a = a, b = b, TE = TE, NIE = a * b, g = (a * b) / TE)
}

## Leading amplitude C_0 in g_m = g* + C_0/m + O(m^-2), and the fixed point g*.
rg_C0 <- function(p) {
  tb <- p$ab * p$bb + p$cb
  tw_val <- p$aw * p$bw + p$cw
  gb <- p$ab * p$bb / tb
  A_star <- (p$aw - p$ab) * p$tw / p$tb
  T_star <- (tw_val - tb) * p$tw / p$tb
  da <- p$aw - p$ab
  db_ <- p$bw - p$bb
  dc <- p$cw - p$cb
  B_star <- (db_ * p$sMw + p$tw * da * (p$aw * db_ + dc)) / p$sMb
  list(gb = gb,
       C0 = ((p$ab * B_star + p$bb * A_star) * tb - p$ab * p$bb * T_star) / tb^2)
}

## Simulate one two-level sample: n_cluster clusters of m_max units each.
rg_simulate <- function(p, n_cluster, m_max, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- n_cluster * m_max
  cl <- rep(seq_len(n_cluster), each = m_max)
  Ab <- stats::rnorm(n_cluster, sd = sqrt(p$tb))[cl]
  Aw <- stats::rnorm(n, sd = sqrt(p$tw))
  Mb <- stats::rnorm(n_cluster, sd = sqrt(p$sMb))[cl]
  Mw <- stats::rnorm(n, sd = sqrt(p$sMw))
  Yb <- stats::rnorm(n_cluster, sd = sqrt(p$sYb))[cl]
  Yw <- stats::rnorm(n, sd = sqrt(p$sYw))
  A <- Ab + Aw
  M <- p$ab * Ab + p$aw * Aw + Mb + Mw
  Y <- p$bb * (p$ab * Ab + Mb) + p$bw * (p$aw * Aw + Mw) +
    p$cb * Ab + p$cw * Aw + Yb + Yw
  data.frame(A = A, M = M, Y = Y, cluster = cl)
}

## The manuscript's canonical flowing and fixed-line parameter sets.
RG_DGP <- list(
  heterogeneous = list(tb = 1, tw = 1, ab = 0.8, aw = 0.5, bb = 0.8, bw = 0.5,
                       cb = 0.1, cw = 0.4, sMb = 0.5, sMw = 0.5,
                       sYb = 0.5, sYw = 0.5),
  homogeneous   = list(tb = 1, tw = 1, ab = 0.6, aw = 0.6, bb = 0.6, bw = 0.6,
                       cb = 0.3, cw = 0.3, sMb = 0.5, sMw = 0.7,
                       sYb = 0.5, sYw = 0.9)
)
