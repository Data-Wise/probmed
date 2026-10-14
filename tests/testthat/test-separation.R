# Separation detector (repair task T11): deterministic LP, not glm warnings.

test_that("complete and quasi-complete separation are detected", {
  skip_if_not_installed("lpSolve")
  x <- c(-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2)
  X <- cbind(1, x)
  expect_identical(.detect_separation(X, as.numeric(x > 0)), "separated")
  # Quasi: the two classes touch at x = 0 only.
  xq <- c(-2, -1, 0, 0, 1, 2)
  expect_identical(.detect_separation(cbind(1, xq), c(0, 0, 0, 1, 1, 1)), "separated")
  # A single class is separated trivially.
  expect_identical(.detect_separation(X, rep(1, 8)), "separated")
})

test_that("overlapping data are not flagged, even with extreme fitted probabilities", {
  skip_if_not_installed("lpSolve")
  set.seed(4)
  x <- stats::rnorm(400)
  y <- stats::rbinom(400, 1, stats::plogis(-0.2 + 0.8 * x))
  expect_identical(.detect_separation(cbind(1, x), y), "ok")
  # Steep but overlapping: fitted probabilities reach ~0/1, the classes still
  # overlap, so this is not separation. One flipped label in each tail.
  xs <- seq(-6, 6, length.out = 200)
  ys <- as.numeric(xs > 0)
  ys[c(5, 196)] <- 1 - ys[c(5, 196)]
  expect_identical(.detect_separation(cbind(1, xs), ys), "ok")
})

test_that("a separated original fit stops pmed() with an informative error", {
  skip_if_not_installed("lpSolve")
  set.seed(5)
  n <- 120
  X <- stats::rbinom(n, 1, 0.5)
  M <- 0.6 * X + stats::rnorm(n)
  d <- data.frame(X = X, M = M, Y = as.numeric(M > 0.3))
  expect_error(suppressWarnings(
    pmed(Y ~ X + M, formula_m = M ~ X, data = d, treatment = "X",
         mediator = "M", family_y = stats::binomial(), method = "plugin")),
    "separation in the outcome model"
  )
})

test_that("non-separated binary fits record their separation status", {
  skip_if_not_installed("lpSolve")
  set.seed(6)
  n <- 200
  X <- stats::rbinom(n, 1, 0.5)
  M <- 0.6 * X + stats::rnorm(n)
  d <- data.frame(X = X, M = M, Y = stats::rbinom(n, 1, stats::plogis(0.8 * M)))
  r <- pmed(Y ~ X + M, formula_m = M ~ X, data = d, treatment = "X",
            mediator = "M", family_y = stats::binomial(),
            method = "nonparametric_bootstrap", n_boot = 10, n_sim = 2000, seed = 1)
  expect_identical(r@diagnostics$separation, c(outcome = "ok"))
  expect_identical(r@diagnostics$separated_draws, 0L)
})

test_that("the detector degrades to 'unchecked' without lpSolve", {
  testthat::local_mocked_bindings(
    requireNamespace = function(...) FALSE, .package = "base"
  )
  expect_identical(.detect_separation(cbind(1, 1:4), c(0, 0, 1, 1)), "unchecked")
})
