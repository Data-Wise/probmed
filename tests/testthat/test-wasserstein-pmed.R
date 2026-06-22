## test-wasserstein-pmed.R — production tests for P_med^W (Paper 4).
## Mirrors the verified facts in pmed-modern/04-wasserstein-pmed/sims/wasserstein_verify.R.

## FACT1-2: w2_1d recovers Gaussian W2 (tested implicitly via pmedW_engine)

test_that("potIF matches the Gaussian closed-form 1-D W2^2 influence function", {
  set.seed(1)
  m <- 0.5; s <- 1.5
  ya <- rnorm(2e5, 0, 1); yb <- rnorm(2e5, m, s)
  if_num <- probmed:::.potIF(ya, yb)
  if_cf  <- (1 - s) * ya^2 - 2 * m * ya - (1 - s)
  if_cf  <- if_cf - mean(if_cf)
  expect_lt(mean(abs(if_num - if_cf)), 0.05)
})

test_that("kill-test: pure variance mediation gives P_med^W>0 (mean-scale blind)", {
  set.seed(1)
  mu <- 0.5
  nu11 <- rnorm(2e5, mu, 2.0); nu10 <- rnorm(2e5, mu, 1.2); nu00 <- rnorm(2e5, mu, 1.0)
  e <- probmed:::.pmedW_engine(nu11, nu10, nu00, eif = FALSE)
  # |sd(nu11)-sd(nu10)| / (|sd(nu11)-sd(nu10)| + |sd(nu10)-sd(nu00)|) = 0.8/1.0 = 0.8
  expect_equal(e$pmedW, 0.8, tolerance = 0.02)
  expect_equal(mean(nu11) - mean(nu10), 0, tolerance = 0.02)  # mean NIE ~ 0
})

test_that("reduction: linear-Gaussian DGP P_med^W ~ classical NIE/(NIE+NDE)", {
  skip_on_cran()
  set.seed(2); n <- 4000
  C <- rnorm(n); A <- rbinom(n, 1, plogis(0.3 * C))
  M <- 0.4 + 0.7 * A + 0.5 * C + rnorm(n)
  Y <- 0.2 + 0.6 * A + 0.8 * M + 0.3 * C + rnorm(n)
  r <- wasserstein_pmed(data.frame(A, M, Y, C), covars = "C", n_mc = 20000L, method = "eif")
  classical <- (0.7 * 0.8) / (0.7 * 0.8 + 0.6)  # a*b / (a*b + c') = 0.483
  expect_lt(abs(r@pmedW - classical), 0.05)
  expect_lt(r@synergy, 0.03)  # no A x M -> synergy ~ 0
  expect_true(inherits(r, "probmed::WassersteinPmedResult"))
})

test_that("EIF SE tracks empirical SD (root-n) on Gaussian corner laws", {
  skip_on_cran()
  set.seed(3)
  mC <- c(1.2, 0.4, 0.0); sC <- c(1.8, 1.1, 1.0)
  eif_se_n3e5 <- probmed:::.pmedW_engine(
    rnorm(3e5, mC[1], sC[1]), rnorm(3e5, mC[2], sC[2]),
    rnorm(3e5, mC[3], sC[3]), eif = TRUE)$se
  eif_se_n2000 <- eif_se_n3e5 * sqrt(3e5 / 2000)  # rescale to n=2000
  emp <- sd(replicate(400, {
    probmed:::.pmedW_engine(
      rnorm(2000, mC[1], sC[1]), rnorm(2000, mC[2], sC[2]),
      rnorm(2000, mC[3], sC[3]), eif = FALSE)$pmedW
  }))
  expect_equal(eif_se_n2000 / emp, 1, tolerance = 0.12)
})

test_that("boundary endpoints: identical corners give P_med^W in {0,1}", {
  set.seed(4); base <- rnorm(1e5)
  e0 <- probmed:::.pmedW_engine(base, base, base + 1, eif = FALSE)
  e1 <- probmed:::.pmedW_engine(base + 1, base, base, eif = FALSE)
  expect_equal(e0$pmedW, 0, tolerance = 0.02)  # nu11 == nu10 -> NIE^W = 0
  expect_equal(e1$pmedW, 1, tolerance = 0.02)  # nu10 == nu00 -> NDE^W = 0
})

test_that("bootstrap method: returns valid WassersteinPmedResult with SE", {
  skip_on_cran()
  set.seed(7); n <- 300
  C <- rnorm(n); A <- rbinom(n, 1, 0.5)
  M <- 0.3 * A + 0.5 * C + rnorm(n)
  Y <- 0.3 * A + 0.8 * M + 0.2 * C + rnorm(n)
  r <- wasserstein_pmed(data.frame(A, M, Y, C), covars = "C",
                        method = "bootstrap", n_boot = 50L, seed = 7L)
  expect_true(inherits(r, "probmed::WassersteinPmedResult"))
  expect_equal(r@method, "bootstrap")
  expect_true(!is.na(r@pmedW_se) && r@pmedW_se > 0)
  expect_length(r@pmedW_ci, 2L)
})

test_that("sinkhorn_div recovers d=2 W2^2 for isotropic Gaussians", {
  skip_on_cran()
  set.seed(11)
  mu1 <- c(1.2, 0.8); mu2 <- c(0.4, 0.3); s1 <- 1.6; s2 <- 1.1
  X <- cbind(rnorm(400, mu1[1], s1), rnorm(400, mu1[2], s1))
  Y <- cbind(rnorm(400, mu2[1], s2), rnorm(400, mu2[2], s2))
  truth <- sum((mu1 - mu2)^2) + 2 * (s1 - s2)^2  # closed form for isotropic Gaussian W2^2
  # tight eps=0.05 + large corners: Sinkhorn convergence is approximate at the default iters;
  # this test asserts the *estimate* (below), so the convergence warning is expected and
  # suppressed here -- it stays live for general use.
  got <- suppressWarnings(sinkhorn_div(X, Y, eps = 0.05))
  expect_equal(got, truth, tolerance = 0.3 * truth + 0.1)
})

test_that("pmedW_md: d=2 corner-law P_med^W recovers truth", {
  skip_on_cran()
  set.seed(11)
  mu <- list(c(1.2, 0.8), c(0.4, 0.3), c(0, 0)); sg <- c(1.6, 1.1, 1.0)
  drw <- function(k, N = 350) cbind(rnorm(N, mu[[k]][1], sg[k]), rnorm(N, mu[[k]][2], sg[k]))
  w2iso <- function(m1, s1, m2, s2) sqrt(sum((m1 - m2)^2) + 2 * (s1 - s2)^2)
  truth <- w2iso(mu[[1]], sg[1], mu[[2]], sg[2]) /
           (w2iso(mu[[1]], sg[1], mu[[2]], sg[2]) + w2iso(mu[[2]], sg[2], mu[[3]], sg[3]))
  r <- suppressWarnings(pmedW_md(drw(1), drw(2), drw(3), eps = 0.1))  # approx Sinkhorn at large n
  expect_equal(r$pmedW, truth, tolerance = 0.04)
})

test_that("wasserstein_pmed(method='dr') runs the user-facing DR path", {
  skip_on_cran()
  set.seed(5); n <- 1400
  C <- rnorm(n); A <- rbinom(n, 1, plogis(0.4 * C))
  M <- 0.3 + 0.6 * A + 0.5 * C + rnorm(n)
  Y <- 0.2 + 0.5 * A + 0.7 * M + 0.3 * C + rnorm(n)
  r <- wasserstein_pmed(data.frame(A, M, Y, C), covars = "C", method = "dr",
                        K = 5L, G = 120L, n_boot = 15L, seed = 5L)
  expect_true(inherits(r, "probmed::WassersteinPmedResult"))
  expect_equal(r@method, "dr")
  expect_true(is.finite(r@pmedW) && r@pmedW >= 0 && r@pmedW <= 1)
  expect_true(!is.na(r@pmedW_se) && r@pmedW_se >= 0)
  expect_length(r@pmedW_ci, 2L)
})

test_that("wasserstein_pmed handles a binary outcome (binY=TRUE)", {
  skip_on_cran()
  set.seed(8); n <- 1200
  C <- rnorm(n); A <- rbinom(n, 1, plogis(0.3 * C))
  M <- 0.4 + 0.6 * A + 0.5 * C + rnorm(n)
  Y <- rbinom(n, 1, plogis(-0.2 + 0.5 * A + 0.7 * M + 0.3 * C))
  r <- wasserstein_pmed(data.frame(A, M, Y, C), covars = "C", binY = TRUE,
                        method = "bootstrap", n_boot = 30L, n_mc = 5000L, seed = 8L)
  expect_true(inherits(r, "probmed::WassersteinPmedResult"))
  expect_true(is.finite(r@pmedW) && r@pmedW >= 0 && r@pmedW <= 1)
})

test_that("eif method falls back to bootstrap at the NIE^W=0 boundary", {
  skip_on_cran()
  # No A -> M effect => true NIE^W = 0; a strong direct effect makes NDE^W (hence D)
  # large, so the estimated NIE^W/D ratio sits well under the 0.05 boundary threshold.
  # n_mc = 20000 keeps the finite-MC upward bias of the empirical W2 small.
  set.seed(9); n <- 1500
  C <- rnorm(n); A <- rbinom(n, 1, 0.5)
  M <- 0.5 * C + rnorm(n)                     # no A -> M effect => NIE^W ~ 0
  Y <- 1.5 * A + 0.8 * M + 0.3 * C + rnorm(n) # strong direct effect => large NDE^W
  expect_message(
    r <- wasserstein_pmed(data.frame(A, M, Y, C), covars = "C", method = "eif",
                          n_boot = 20L, n_mc = 20000L, seed = 9L),
    "boundary"
  )
  expect_true(grepl("boundary", r@method))
  expect_true(r@boundary)
})

test_that("validator warns when NIE^W is at the non-regular boundary", {
  expect_warning(
    WassersteinPmedResult(
      pmedW = 0, pmedW_ci = c(0, 0), pmedW_se = NA_real_,
      NIE_W = 1e-5, NDE_W = 0.5, TE_W = 0.5, synergy = 0,
      method = "eif", n = 100L, ci_level = 0.95),
    "non-regular boundary"
  )
})

test_that("print method covers boundary and synergy warning branches", {
  r <- WassersteinPmedResult(
    pmedW = 0.5, pmedW_ci = c(0.3, 0.7), pmedW_se = 0.1,
    NIE_W = 0.4, NDE_W = 0.4, TE_W = 0.5, synergy = 0.3,  # > 10% of (NIE+NDE)=0.08
    boundary = TRUE, method = "bootstrap-boundary", n = 100L, ci_level = 0.95)
  out <- capture.output(print(r))
  expect_true(any(grepl("Wasserstein P_med", out)))
  expect_true(any(grepl("boundary", out)))
  expect_true(any(grepl("triangle gap", out)))
  expect_invisible(print(r))
})

test_that("sinkhorn early-convergence break is exercised (large eps)", {
  set.seed(13)
  X <- matrix(rnorm(40), 20, 2); Y <- matrix(rnorm(40) + 5, 20, 2)
  expect_gt(sinkhorn_div(X, Y, eps = 1.0, iters = 400L), 0)
})

test_that("DR corner-law estimator is multiply robust", {
  skip_on_cran()
  set.seed(5); n <- 4000
  C <- rnorm(n); A <- rbinom(n, 1, plogis(0.5 * C))
  M <- 0.3 + 0.6 * A + 0.5 * C + rnorm(n)
  Y <- 0.2 + 0.5 * A + 0.7 * M + 0.4 * C + 0.3 * A * M + rnorm(n)
  d <- data.frame(A, M, Y, C)
  # oracle truth via large MC g-computation
  Ct <- rnorm(2e5)
  cl_truth <- function(a, ap) {
    Mt <- 0.3 + 0.6 * ap + 0.5 * Ct + rnorm(2e5)
    0.2 + 0.5 * a + 0.7 * Mt + 0.4 * Ct + 0.3 * a * Mt + rnorm(2e5)
  }
  U <- (seq_len(2048) - 0.5) / 2048
  q <- function(x) quantile(x, U, names = FALSE)
  N_true <- sqrt(mean((q(cl_truth(1,1)) - q(cl_truth(1,0)))^2))
  D_true <- sqrt(mean((q(cl_truth(1,0)) - q(cl_truth(0,0)))^2))
  truth  <- N_true / (N_true + D_true)
  for (ms in c("none", "outcome", "propensity")) {
    e <- pmedW_dr(d, covars = "C", K = 5L, misspec = ms, seed = 5L)
    expect_lt(abs(e$pmedW - truth), 0.05)
  }
})

test_that("DR estimator is robust to outcome- or propensity-model misspecification", {
  skip_on_cran()
  set.seed(11); n <- 1500
  C <- rnorm(n); A <- rbinom(n, 1, plogis(0.3 * C))
  M <- 0.4 + 0.7 * A + 0.5 * C + rnorm(n)
  Y <- 0.2 + 0.6 * A + 0.8 * M + 0.3 * C + rnorm(n)
  d <- data.frame(A, M, Y, C)
  truth <- (0.7 * 0.8) / (0.7 * 0.8 + 0.6)             # classical P_med ~ 0.483
  none <- probmed:::pmedW_dr(d, covars = "C", K = 5L, misspec = "none",       G = 100L, seed = 1L)$pmedW
  mo   <- probmed:::pmedW_dr(d, covars = "C", K = 5L, misspec = "outcome",    G = 100L, seed = 1L)$pmedW
  mp   <- probmed:::pmedW_dr(d, covars = "C", K = 5L, misspec = "propensity", G = 100L, seed = 1L)$pmedW
  expect_lt(abs(none - truth), 0.08)
  ## multiple robustness: stays near truth when ONE nuisance is wrong
  expect_lt(abs(mo - truth), 0.12)
  expect_lt(abs(mp - truth), 0.12)
})

test_that("near the NIE^W = 0 boundary, eif switches to a bootstrap CI and flags it", {
  skip_on_cran()
  set.seed(12); n <- 400
  C <- rnorm(n); A <- rbinom(n, 1, 0.5)
  M <- 0.3 + 0.8 * A + 0.4 * C + rnorm(n)
  Y <- 0.5 * A + 0.3 * C + rnorm(n)               # M absent from Y => no mediation => NIE^W ~ 0
  expect_message(
    r <- wasserstein_pmed(data.frame(A, M, Y, C), covars = "C",
                          method = "eif", n_boot = 50L, seed = 12L),
    "boundary")
  expect_true(r@boundary)
  expect_true(!is.na(r@pmedW_se) && r@pmedW_se > 0)   # bootstrap CI still produced
  expect_length(r@pmedW_ci, 2L)
})

test_that("sinkhorn warns when it fails to converge (too few iterations)", {
  set.seed(13)
  X <- matrix(rnorm(60), ncol = 2); Y <- matrix(rnorm(60) + 3, ncol = 2)
  expect_warning(probmed:::.sinkhorn_cost(X, Y, eps = 0.01, iters = 1L), "not converged")
})
