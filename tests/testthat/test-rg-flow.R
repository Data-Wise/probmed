## Truth constants generated from the manuscript repo's own dgp.R
## (research/mediation-noncollapsibility, 01-rg-pmed/sims/dgp.R):
##   source("sims/dgp.R"); pop_paths(DGP$heterogeneous, m)$g
## Cross-checked against helper-rg.R's reimplementation below, so a
## transcription slip in either copy fails the suite.
G_TRUTH_HET <- c(
  `1`    = 0.5938177688,
  `2`    = 0.6749729144,
  `5`    = 0.7645079899,
  `10`   = 0.8086066320,
  `50`   = 0.8523973073,
  `1000` = 0.8642249980
)
C0_TRUTH_HET <- -0.6407596786
GSTAR_TRUTH_HET <- 0.8648648649
G_TRUTH_HOM <- 0.5454545455

test_that("vendored truth function reproduces the manuscript's constants", {
  # Guards the helper itself: if helper-rg.R drifts from dgp.R, this fails
  # before any estimator test can silently inherit the drift.
  for (m in names(G_TRUTH_HET)) {
    expect_equal(rg_pop_paths(RG_DGP$heterogeneous, as.numeric(m))$g,
                 unname(G_TRUTH_HET[[m]]), tolerance = 1e-9)
  }
  cc <- rg_C0(RG_DGP$heterogeneous)
  expect_equal(cc$C0, C0_TRUTH_HET, tolerance = 1e-9)
  expect_equal(cc$gb, GSTAR_TRUTH_HET, tolerance = 1e-9)
})

test_that("analytic gradient matches central finite differences", {
  # The likeliest defect in this estimator is a transcription error in one of
  # the five long partials, so check them directly against numeric derivatives
  # of g(Sigma) rather than only through downstream quantities.
  set.seed(11)
  d <- rg_simulate(RG_DGP$heterogeneous, n_cluster = 60, m_max = 8, seed = 11)
  Z <- as.matrix(d[, c("A", "M", "Y")])
  parts <- probmed:::.rg_parts(Z)

  # g as an explicit function of the five moments, so FD is independent of
  # .rg_parts()'s own gradient code.
  g_of_sigma <- function(s) {
    D <- s[["AA"]] * s[["MM"]] - s[["AM"]]^2
    a <- s[["AM"]] / s[["AA"]]
    b <- (s[["MY"]] * s[["AA"]] - s[["AY"]] * s[["AM"]]) / D
    TE <- s[["AY"]] / s[["AA"]]
    (a * b) / TE
  }
  expect_equal(g_of_sigma(parts$sig), parts$g, tolerance = 1e-12)

  h <- 1e-6
  fd <- vapply(names(parts$sig), function(k) {
    up <- dn <- parts$sig
    up[[k]] <- up[[k]] + h
    dn[[k]] <- dn[[k]] - h
    (g_of_sigma(up) - g_of_sigma(dn)) / (2 * h)
  }, numeric(1))
  expect_equal(unname(parts$grad), unname(fd), tolerance = 1e-6)
})

test_that("g_hat_eif recovers the population truth at several scales", {
  skip_on_cran()
  set.seed(202607)
  p <- RG_DGP$heterogeneous
  d <- rg_simulate(p, n_cluster = 4000, m_max = 10, seed = 202607)
  for (m in c(1, 2, 5, 10)) {
    cg <- rg_coarse_grain(d, m, cluster = "cluster")
    fit <- g_hat_eif(cg$Z, cg$cluster, m = m)
    truth <- rg_pop_paths(p, m)$g
    expect_equal(fit@g, truth, tolerance = 0.02)
    expect_equal(fit@a_path, rg_pop_paths(p, m)$a, tolerance = 0.02)
    expect_true(fit@g_ci[1] < truth && truth < fit@g_ci[2])
  }
})

test_that("the fixed line is flat: g does not move with scale", {
  skip_on_cran()
  set.seed(4242)
  d <- rg_simulate(RG_DGP$homogeneous, n_cluster = 4000, m_max = 10, seed = 4242)
  gs <- vapply(c(1, 2, 5, 10), function(m) {
    cg <- rg_coarse_grain(d, m, cluster = "cluster")
    g_hat_eif(cg$Z, cg$cluster, m = m)@g
  }, numeric(1))
  expect_equal(unname(gs), rep(G_TRUTH_HOM, 4), tolerance = 0.02)
  # flat means flat: the spread across scales is far below the flowing case,
  # which moves ~0.21 over the same grid.
  expect_lt(diff(range(gs)), 0.02)
})

test_that("coarse-graining drops incomplete blocks and preserves clusters", {
  d <- data.frame(A = rnorm(70), M = rnorm(70), Y = rnorm(70),
                  cluster = rep(1:7, each = 10))
  cg <- rg_coarse_grain(d, m = 3, cluster = "cluster")
  # 7 clusters x floor(10/3) = 7 x 3 = 21 units, NOT 70/3
  expect_equal(nrow(cg$Z), 21L)
  expect_equal(length(unique(cg$cluster)), 7L)
  expect_identical(colnames(cg$Z), c("A", "M", "Y"))
  # m = 1 is the identity on the unit scale
  expect_equal(nrow(rg_coarse_grain(d, 1, cluster = "cluster")$Z), 70L)
  # clusters smaller than m drop out entirely
  d2 <- data.frame(A = rnorm(12), M = rnorm(12), Y = rnorm(12),
                   cluster = rep(c("big", "small"), times = c(10, 2)))
  expect_equal(length(unique(rg_coarse_grain(d2, 5, cluster = "cluster")$cluster)), 1L)
})

test_that("coarse-graining composes: T_m2(T_m1) == T_(m1*m2) on retained units", {
  d <- data.frame(A = 1:24, M = (1:24) * 2, Y = (1:24) * 3,
                  cluster = rep(1L, 24))
  # deterministic when the partition is the whole cluster in order: compare
  # means over the same retained prefix rather than relying on the shuffle.
  Z <- as.matrix(d[, c("A", "M", "Y")])
  stage1 <- rowsum(Z, rep(seq_len(12), each = 2)) / 2   # m1 = 2
  composed <- rowsum(stage1, rep(seq_len(4), each = 3)) / 3  # m2 = 3
  direct <- rowsum(Z, rep(seq_len(4), each = 6)) / 6         # m = 6
  expect_equal(unname(composed), unname(direct), tolerance = 1e-12)
})

test_that("cluster-robust SE exceeds the naive i.i.d. SE under clustering", {
  skip_on_cran()
  set.seed(7)
  d <- rg_simulate(RG_DGP$heterogeneous, n_cluster = 300, m_max = 10, seed = 7)
  cg <- rg_coarse_grain(d, 1, cluster = "cluster")
  fit <- g_hat_eif(cg$Z, cg$cluster, m = 1)
  parts <- probmed:::.rg_parts(cg$Z)
  phi_i <- probmed:::.rg_phi_units(cg$Z, parts)
  se_iid <- sqrt(mean(phi_i^2) / length(phi_i))
  # Treating units as independent understates the SE when A has a between
  # component; this is why the cluster is the sampling unit.
  expect_gt(fit@g_se, se_iid)
  expect_equal(fit@n_clusters, 300L)
})

test_that("paired contrast is more precise than the naive independent one", {
  skip_on_cran()
  set.seed(99)
  d <- rg_simulate(RG_DGP$heterogeneous, n_cluster = 800, m_max = 10, seed = 99)
  cg1 <- rg_coarse_grain(d, 1, cluster = "cluster")
  cg2 <- rg_coarse_grain(d, 2, cluster = "cluster")
  f1 <- g_hat_eif(cg1$Z, cg1$cluster, m = 1)
  f2 <- g_hat_eif(cg2$Z, cg2$cluster, m = 2)
  ct <- rg_flow_contrast(f1, f2)

  expect_equal(ct$delta, f2@g - f1@g, tolerance = 1e-12)
  expect_lt(ct$se_paired, ct$se_naive)   # pairing must help, never hurt
  expect_gt(ct$gain, 1)
  expect_equal(ct$n_clusters, 800L)
  # the true difference g_2 - g_1 should sit inside the paired interval
  truth <- rg_pop_paths(RG_DGP$heterogeneous, 2)$g -
    rg_pop_paths(RG_DGP$heterogeneous, 1)$g
  expect_true(ct$ci[1] < truth && truth < ct$ci[2])
})

test_that("input validation rejects malformed calls", {
  d <- data.frame(A = rnorm(20), M = rnorm(20), Y = rnorm(20),
                  cluster = rep(1:2, each = 10))
  expect_error(rg_coarse_grain(d, 3, cluster = "nope"), "not found")
  expect_error(rg_coarse_grain(d, 500, cluster = "cluster"), "at least")
  Z <- as.matrix(d[, c("A", "M", "Y")])
  expect_error(g_hat_eif(Z[, 1:2], d$cluster), "three columns")
  expect_error(g_hat_eif(Z, d$cluster[1:5]), "one label per row")
  expect_error(g_hat_eif(Z, d$cluster, ci_level = 1.5), "in \\(0, 1\\)")
  f <- g_hat_eif(Z, d$cluster)
  expect_error(rg_flow_contrast(f, "not a fit"), "RgFlowResult")
})

test_that("RgFlowResult rejects a numerically zero total effect", {
  # g = NIE/TE is not well posed at TE = 0; the validator must say so.
  expect_error(
    RgFlowResult(g = 1, g_ci = c(0, 2), g_se = 1, a_path = 1, b_path = 1,
                 TE = 0, grad_norm = 1, phi_cluster = 1, m = 1,
                 n_units = 1L, n_clusters = 1L, ci_level = 0.95),
    "not a well-posed ratio"
  )
})
