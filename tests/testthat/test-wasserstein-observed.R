test_that("public Wasserstein inference uses observed rows, not simulated draws", {
  set.seed(20261020)
  n <- 500
  C <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.3 * C))
  M <- 0.7 * A + 0.5 * C + rnorm(n)
  Y <- 0.6 * A + 0.8 * M + 0.3 * C + rnorm(n)
  d <- data.frame(A, M, Y, C)
  a <- wasserstein_pmed(d, n_mc = 2000L, seed = 17L)
  b <- wasserstein_pmed(d, n_mc = 200000L, seed = 17L)
  expect_identical(a@pmedW, b@pmedW)
  expect_identical(a@pmedW_se, b@pmedW_se)
  expect_gt(a@pmedW_se, 0.02)
})

wasserstein_fixture <- function(binary = FALSE, n = 400L) {
  set.seed(810 + binary)
  C <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.3 * C))
  M <- 0.7 * A + 0.5 * C + rnorm(n)
  eta <- -0.5 + 0.6 * A + 0.8 * M + 0.2 * A * M + 0.3 * C
  Y <- if (binary) rbinom(n, 1, plogis(eta)) else eta + rnorm(n)
  data.frame(A, M, Y, C)
}

test_that("observed-row IF agrees with weighted-refit contamination derivatives", {
  for (binary in c(FALSE, TRUE)) {
    d <- wasserstein_fixture(binary)
    fit <- probmed:::.wasserstein_observed(d, binY = binary)
    n <- nrow(d)
    eps <- 1e-5
    for (i in c(1L, 23L, which.max(d$C), which.min(d$Y))) {
      w <- rep((1 - eps) / n, n)
      w[i] <- w[i] + eps
      # Refits, including weighted residual variances and empirical C, rather
      # than perturbing a fitted-model coefficient in isolation.
      shifted <- suppressWarnings(probmed:::.wasserstein_observed(
        d,
        binY = binary, influence = FALSE, weights = w
      ))
      numerical <- (shifted$pmedW - fit$pmedW) / eps
      expect_equal(unname(fit$influence[i]), numerical, tolerance = 0.005)
    }
    expect_equal(mean(fit$influence), 0, tolerance = 1e-7)
    refined <- probmed:::.wasserstein_observed(d, binY = binary, order = 256L)
    expect_lt(abs(fit$pmedW - refined$pmedW), 1e-4)
    expect_lt(abs(fit$se / refined$se - 1), 0.01)
  }
})

test_that("observed-data sandwich agrees with row bootstrap on regular fixtures", {
  skip_on_cran()
  for (binary in c(FALSE, TRUE)) {
    d <- wasserstein_fixture(binary)
    r <- wasserstein_pmed(d, binY = binary)
    # The initial B=400 binary comparison differed by 22.5%; increasing B to
    # 2,000 (same fixture, seed and unchanged tolerance) gave 18.1%. Retain the
    # higher-precision check rather than relaxing the prespecified tolerance.
    bs <- wasserstein_pmed(d,
      binY = binary, method = "bootstrap",
      n_boot = if (binary) 2000L else 400L
    )
    expect_lt(abs(r@pmedW_se / bs@pmedW_se - 1), 0.20)
    expect_equal(r@pmedW, bs@pmedW)
  }
})

test_that("Gaussian mixture transport agrees with an independent CDF integration", {
  d <- wasserstein_fixture(n = 200L)
  fm <- lm(M ~ A + C, d)
  fy <- lm(Y ~ A * M + C, d)
  vm <- mean(resid(fm)^2)
  vy <- mean(resid(fy)^2)
  # Independent stats::integrate over u, rather than production Gaussian rules.
  quant <- function(u, a, ap) {
    mm <- predict(fm, transform(d, A = ap))
    intercept <- predict(fy, transform(d, A = a, M = 0))
    slope <- predict(fy, transform(d, A = a, M = 1)) - intercept
    mu <- intercept + slope * mm
    s <- sqrt(vy + slope^2 * vm)
    vapply(u, function(p) {
      ends <- range(mu + s * qnorm(p))
      uniroot(function(y) mean(pnorm(y, mu, s)) - p, ends, tol = 1e-10)$root
    }, numeric(1))
  }
  distance <- function(a, ap, b, bp) {
    sqrt(integrate(
      function(u) (quant(u, a, ap) - quant(u, b, bp))^2,
      0, 1,
      rel.tol = 1e-6
    )$value)
  }
  nie <- distance(1, 1, 1, 0)
  nde <- distance(1, 0, 0, 0)
  fit <- probmed:::.wasserstein_observed(d)
  expect_equal(fit$pmedW, nie / (nie + nde), tolerance = 1e-4)
})

test_that("Bernoulli corner transport uses absolute probability differences", {
  d <- wasserstein_fixture(TRUE)
  fm <- lm(M ~ A + C, d)
  fy <- glm(Y ~ A * M + C, d, family = binomial())
  vm <- mean(resid(fm)^2)
  p <- function(a, ap) {
    mu <- predict(fm, transform(d, A = ap))
    base <- predict(fy, transform(d, A = a, M = 0))
    slope <- predict(fy, transform(d, A = a, M = 1)) - base
    mean(vapply(seq_len(nrow(d)), function(i) {
      integrate(
        function(z) plogis(base[i] + slope[i] * (mu[i] + sqrt(vm) * z)) * dnorm(z),
        -Inf, Inf,
        rel.tol = 1e-9
      )$value
    }, numeric(1)))
  }
  ps <- c(p(1, 1), p(1, 0), p(0, 0))
  truth <- sqrt(abs(ps[1] - ps[2])) /
    (sqrt(abs(ps[1] - ps[2])) + sqrt(abs(ps[2] - ps[3])))
  expect_equal(probmed:::.wasserstein_observed(d, binY = TRUE)$pmedW,
    truth,
    tolerance = 1e-7
  )
})

test_that("observed-data inference rejects invalid or unidentifiable fits", {
  d <- wasserstein_fixture()
  bad <- d
  bad$M[1] <- NA_real_
  expect_error(probmed:::.wasserstein_observed(bad), "complete finite")
  bad <- d
  bad$A <- 0
  expect_error(probmed:::.wasserstein_observed(bad), "binary A")
  expect_error(probmed:::.wasserstein_observed(d, binY = TRUE), "binary Y")
  bad <- d
  bad$C <- bad$A
  expect_error(probmed:::.wasserstein_observed(bad), "full rank")
  bad <- d
  bad$Y <- bad$A + bad$M + bad$C
  expect_error(probmed:::.wasserstein_observed(bad), "positive residual")
  expect_error(probmed:::.wasserstein_observed(d, weights = rep(1, nrow(d))), "weights")
  w <- rep(1 / nrow(d), nrow(d))
  w[1:2] <- w[1:2] + c(0.001, -0.001)
  expect_error(probmed:::.wasserstein_observed(d, weights = w), "equal observed")
})

test_that("both zero-component boundaries are treated as nonregular", {
  d <- wasserstein_fixture(n = 200L)
  fy <- coef(lm(Y ~ A * M + C, d))
  d$Y <- d$Y - fy["A"] * d$A - fy["A:M"] * d$A * d$M
  expect_message(r <- wasserstein_pmed(d, n_boot = 10L), "boundary")
  expect_true(r@boundary)
  expect_equal(r@pmedW, 1, tolerance = 1e-7)
  expect_identical(r@method, "bootstrap-boundary")
  fm <- coef(lm(M ~ A + C, d))
  d$M <- d$M - fm["A"] * d$A
  # Refit Y without A or A:M after changing M: the fitted corners all coincide.
  fy <- coef(lm(Y ~ A * M + C, d))
  d$Y <- d$Y - fy["A"] * d$A - fy["A:M"] * d$A * d$M
  expect_error(probmed:::.wasserstein_observed(d), "both components vanish")
})
