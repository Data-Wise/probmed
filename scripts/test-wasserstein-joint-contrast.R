source("scripts/wasserstein-joint-contrast.R")
testthat::test_that("projection handles axes and undefined origin without a point convention", {
  testthat::expect_equal(wj_interval(c(0, 0), diag(2)), c(0, 1))
  testthat::expect_equal(wj_interval(c(0.1, 0.1), diag(2)), c(0, 1))
  testthat::expect_equal(wj_interval(c(0, 2), diag(0.01, 2))[1], 0)
  testthat::expect_equal(wj_interval(c(2, 0), diag(0.01, 2))[2], 1)
  testthat::expect_error(wj_interval(c(1, 1), matrix(1, 2, 2)), "positive definite")
  testthat::expect_error(wj_interval(c(1, 1), matrix(c(1, 0, 1, 1), 2)), "symmetric")
  testthat::expect_error(wj_interval(c(1, NA), diag(2)), "finite")
})
testthat::test_that("projection agrees with independent dense ellipse boundaries", {
  v <- matrix(c(0.04, 0.012, 0.012, 0.09), 2)
  angle <- seq(0, 2 * pi, length.out = 200001L)
  boundary <- sqrt(stats::qchisq(0.95, 2)) *
    t(chol(v)) %*% rbind(cos(angle), sin(angle))
  for (d in list(c(2, 3), c(-2, 3), c(0.1, 3), c(3, 0.1))) {
    x <- boundary + d
    p <- sqrt(abs(x[1, ])) / (sqrt(abs(x[1, ])) + sqrt(abs(x[2, ])))
    expected <- range(p)
    if (d[1]^2 / v[1, 1] <= stats::qchisq(0.95, 2)) expected[1] <- 0
    if (d[2]^2 / v[2, 2] <= stats::qchisq(0.95, 2)) expected[2] <- 1
    ci <- wj_interval(d, v)
    testthat::expect_equal(ci, expected, tolerance = 1e-7)
    testthat::expect_true(ci[1] <= wb_ratio(d) && ci[2] >= wb_ratio(d))
  }
})
testthat::test_that("projection preserves scale, signs and complement symmetry", {
  d <- c(0.2, 0.3)
  v <- matrix(c(0.001, 0.0002, 0.0002, 0.002), 2)
  ci <- wj_interval(d, v)
  testthat::expect_equal(wj_interval(100 * d, 10000 * v), ci, tolerance = 1e-9)
  signs <- diag(c(-1, 1))
  testthat::expect_equal(wj_interval(drop(signs %*% d), signs %*% v %*% signs), ci)
  testthat::expect_equal(wj_interval(rev(d), v[2:1, 2:1]), 1 - rev(ci), tolerance = 1e-9)
  testthat::expect_true(all(wj_interval(d, v, level = 0.99) * c(-1, 1) >= ci * c(-1, 1)))
})
