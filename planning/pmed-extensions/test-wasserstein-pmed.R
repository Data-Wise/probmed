## test-wasserstein-pmed.R — prototype tests for P_med^W (testthat ed.3).
## Run: Rscript -e 'testthat::test_file("test-wasserstein-pmed.R")'
## Mirrors the verified facts in 04-wasserstein-pmed/sims/wasserstein_verify.R.
library(testthat)
source("wasserstein_pmed.R")

test_that("potIF matches the Gaussian closed-form 1-D W2^2 influence function", {
  set.seed(1)
  m <- 0.5; s <- 1.5
  ya <- rnorm(2e5, 0, 1); yb <- rnorm(2e5, m, s)
  if_num <- .potIF(ya, yb)
  if_cf  <- (1 - s) * ya^2 - 2 * m * ya - (1 - s)
  if_cf  <- if_cf - mean(if_cf)
  expect_lt(mean(abs(if_num - if_cf)), 0.05)          # sampling + grid error
})

test_that("kill-test: pure variance mediation gives P_med^W>0 (mean scale blind)", {
  set.seed(1)
  mu <- 0.5
  nu11 <- rnorm(2e5, mu, 2.0); nu10 <- rnorm(2e5, mu, 1.2); nu00 <- rnorm(2e5, mu, 1.0)
  e <- .pmedW_engine(nu11, nu10, nu00, eif = FALSE)
  expect_equal(e$pmedW, 0.8, tolerance = 0.02)         # |2-1.2|/(|2-1.2|+|1.2-1|)
  expect_equal(mean(nu11) - mean(nu10), 0, tolerance = 0.02)  # mean NIE ~ 0
})

test_that("reduction: additive linear-Gaussian DGP -> P_med^W ~ classical NIE/(NIE+NDE)", {
  set.seed(2); n <- 4000
  C <- rnorm(n); A <- rbinom(n, 1, plogis(0.3 * C))
  M <- 0.4 + 0.7 * A + 0.5 * C + rnorm(n)
  Y <- 0.2 + 0.6 * A + 0.8 * M + 0.3 * C + rnorm(n)
  r <- wasserstein_pmed(data.frame(A, M, Y, C), covars = "C", n_mc = 20000L, method = "eif")
  classical <- (0.7 * 0.8) / (0.7 * 0.8 + 0.6)         # a*b / (a*b + c') = 0.483 (equal-var corners)
  # absolute tolerance: plug-in W2 is positively biased at finite n (the EIF/debiasing + sim
  # study address this); at n=4000 the gap is ~0.02, well within a plug-in's finite-sample error.
  expect_lt(abs(r@pmedW - classical), 0.04)
  expect_lt(r@synergy, 0.03)                           # no A x M => synergy ~ 0
})

test_that("EIF SE tracks empirical SD (root-n) on Gaussian corner laws", {
  set.seed(3)
  mC <- c(1.2, 0.4, 0.0); sC <- c(1.8, 1.1, 1.0)
  eif <- .pmedW_engine(rnorm(3e5, mC[1], sC[1]), rnorm(3e5, mC[2], sC[2]),
                       rnorm(3e5, mC[3], sC[3]), eif = TRUE)$se
  # rescale EIF SE (computed at n=3e5) to n=2000 for comparison
  eif2000 <- eif * sqrt(3e5 / 2000)
  emp <- sd(replicate(400, {
    .pmedW_engine(rnorm(2000, mC[1], sC[1]), rnorm(2000, mC[2], sC[2]),
                  rnorm(2000, mC[3], sC[3]), eif = FALSE)$pmedW
  }))
  expect_equal(eif2000 / emp, 1, tolerance = 0.12)
})

test_that("boundary endpoints: identical corners give P_med^W in {0,1}", {
  set.seed(4); base <- rnorm(1e5)
  e0 <- .pmedW_engine(base, base, base + 1, eif = FALSE)   # nu11==nu10 -> NIE^W=0
  e1 <- .pmedW_engine(base + 1, base, base, eif = FALSE)   # nu10==nu00 -> NDE^W=0
  expect_equal(e0$pmedW, 0, tolerance = 0.02)
  expect_equal(e1$pmedW, 1, tolerance = 0.02)
})

test_that("entropic Sinkhorn divergence recovers d=2 corner-law P_med^W", {
  set.seed(11)
  mu <- list(c(1.2, 0.8), c(0.4, 0.3), c(0, 0)); sg <- c(1.6, 1.1, 1.0)
  drw <- function(k, N = 350) cbind(rnorm(N, mu[[k]][1], sg[k]), rnorm(N, mu[[k]][2], sg[k]))
  w2iso <- function(m1, s1, m2, s2) sqrt(sum((m1 - m2)^2) + 2 * (s1 - s2)^2)
  truth <- w2iso(mu[[1]], sg[1], mu[[2]], sg[2]) /
           (w2iso(mu[[1]], sg[1], mu[[2]], sg[2]) + w2iso(mu[[2]], sg[2], mu[[3]], sg[3]))
  r <- pmedW_md(drw(1), drw(2), drw(3), eps = 0.1)
  expect_equal(r$pmedW, truth, tolerance = 0.04)   # debiased Sinkhorn, finite n/eps
})
