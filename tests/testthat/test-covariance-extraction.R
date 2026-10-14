# Covariance-extraction validation (repair task T8).
#
# medfit's joint vcov is block diagonal with each block equal to vcov(glm).
# That is exact for the sequential likelihood, which factorizes into the
# mediator and outcome models, so the Fisher information has no cross block.
# It is a model-based result, not a claim that the two estimates are
# independent under misspecification. The independent reference is a stacked
# estimating-equation (sandwich) covariance written from scratch, with no
# probmed or medfit helper.

cov_fixture <- function(binary_y = FALSE, shared_cov = TRUE, n = 300, seed = 21) {
  set.seed(seed)
  X <- stats::rbinom(n, 1, 0.5)
  C <- stats::rnorm(n)
  Z <- stats::rnorm(n)
  M <- 0.5 * X + 0.4 * C + stats::rnorm(n)
  eta <- 0.3 * X + 0.4 * M + (if (shared_cov) 0.4 * C else 0.4 * Z)
  Y <- if (binary_y) stats::rbinom(n, 1, stats::plogis(eta)) else eta + stats::rnorm(n)
  list(
    d = data.frame(X = X, M = M, Y = Y, C = C, Z = Z),
    fm = M ~ X + C,
    fy = if (shared_cov) Y ~ X + M + C else Y ~ X + M + Z,
    fam = if (binary_y) stats::binomial() else stats::gaussian()
  )
}

# Independent reference: the stacked estimating-equation (sandwich) covariance of
# (beta_m, beta_y), computed from scratch. Its cross block is O_p(n^-1/2) when
# both equations are correctly specified; it is not zero under misspecification.
sandwich_corr_ab <- function(d, fm, fy, binary_y) {
  Xm <- stats::model.matrix(fm, d)
  Xy <- stats::model.matrix(fy, d)
  m <- d[[all.vars(fm)[1]]]
  y <- d[[all.vars(fy)[1]]]
  em <- m - drop(Xm %*% qr.solve(Xm, m))
  if (binary_y) {
    by <- stats::glm.fit(Xy, y, family = stats::binomial())$coefficients
    mu <- stats::plogis(drop(Xy %*% by))
    ey <- y - mu
    Ay <- crossprod(Xy, Xy * (mu * (1 - mu)))
  } else {
    ey <- y - drop(Xy %*% qr.solve(Xy, y))
    Ay <- crossprod(Xy)
  }
  km <- ncol(Xm)
  psi <- cbind(Xm * em, Xy * ey)
  Ainv <- matrix(0, ncol(psi), ncol(psi))
  Ainv[seq_len(km), seq_len(km)] <- solve(crossprod(Xm))
  Ainv[-seq_len(km), -seq_len(km)] <- solve(Ay)
  V <- Ainv %*% crossprod(psi) %*% Ainv
  ia <- which(colnames(Xm) == "X")
  ib <- km + which(colnames(Xy) == "M")
  V[ia, ib] / sqrt(V[ia, ia] * V[ib, ib])
}

for (cfg in list(list(binary_y = FALSE, shared_cov = TRUE),
                 list(binary_y = TRUE, shared_cov = TRUE),
                 list(binary_y = FALSE, shared_cov = FALSE))) {
  label <- sprintf("binary_y=%s, shared_cov=%s", cfg$binary_y, cfg$shared_cov)
  test_that(paste("extracted vcov matches the factorized likelihood:", label), {
    fx <- do.call(cov_fixture, cfg)
    fit_m <- stats::glm(fx$fm, data = fx$d)
    fit_y <- stats::glm(fx$fy, data = fx$d, family = fx$fam)
    ex <- extract_mediation(fit_m, treatment = "X", mediator = "M",
                            model_y = fit_y, data = fx$d)
    V <- ex@vcov
    mi <- grep("^m_", rownames(V))
    yi <- grep("^y_", rownames(V))
    # Contract: zero cross block, blocks equal to each model's own vcov.
    expect_equal(max(abs(V[mi, yi])), 0)
    expect_equal(unname(V[mi, mi]), unname(stats::vcov(fit_m)), tolerance = 1e-12)
    expect_equal(unname(V[yi, yi]), unname(stats::vcov(fit_y)), tolerance = 1e-12)
    # Aliases duplicate their source coefficients exactly.
    expect_equal(V["a", "a"], V["m_X", "m_X"])
    expect_equal(V["b", "b"], V["y_M", "y_M"])
  })
}

test_that("the sandwich cross-correlation of (a_hat, b_hat) vanishes when correctly specified", {
  # n = 20,000 keeps the O(n^-1/2) sampling noise far below the 0.05 limit
  # (observed |corr| <= 0.016 in these three configurations).
  for (cfg in list(list(binary_y = FALSE, shared_cov = TRUE),
                   list(binary_y = TRUE, shared_cov = TRUE),
                   list(binary_y = FALSE, shared_cov = FALSE))) {
    fx <- do.call(cov_fixture, c(cfg, n = 20000))
    r <- sandwich_corr_ab(fx$d, fx$fm, fx$fy, cfg$binary_y)
    expect_lt(abs(r), 0.05)
  }
})

test_that("the sandwich reference detects a real cross-equation correlation", {
  # Positive control: the outcome depends on (e_m^2 - 1) * X, which the
  # additive outcome model omits. The omitted term correlates the mediator
  # residual with the outcome score, so the zero cross block no longer holds.
  set.seed(8)
  n <- 20000
  X <- stats::rbinom(n, 1, 0.5)
  e <- stats::rnorm(n)
  M <- 0.5 * X + e
  Y <- 0.4 * M + 1.5 * (e^2 - 1) * X + stats::rnorm(n)
  r <- sandwich_corr_ab(data.frame(X = X, M = M, Y = Y), M ~ X, Y ~ X + M,
                        binary_y = FALSE)
  expect_gt(abs(r), 0.2)
})
