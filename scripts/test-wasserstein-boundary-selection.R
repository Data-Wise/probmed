source("scripts/wasserstein-boundary-selection.R")
v <- matrix(c(1, .3, .3, 1), 2)
d <- c(.4, 8)
stopifnot(ws_point(matrix(d, 1), v, qnorm(.975)) == 0)
stopifnot(ws_point(matrix(rev(d), 1), v, qnorm(.975)) == 1)
x <- rbind(c(4, 6), c(-4, 6), c(.4, .3))
stopifnot(max(abs(ws_point(x, v, qnorm(.975)) -
                    apply(x, 1, wb_correct, v = v))) < 1e-14)
set.seed(817)
independent <- matrix(rnorm(262144 * 2), ncol = 2)
set.seed(816)
z <- matrix(rnorm(32768 * 2), ncol = 2)
z <- rbind(z, -z)
for (lambda in qnorm(c(.975, .9975))) {
  for (d in list(c(.4, 8), c(4, 6), c(.4, .3))) {
    a <- ws_fit(d, v, lambda, z)
    b <- ws_fit(d, v, lambda, independent)
    stopifnot(abs(a$se / b$se - 1) < .10)
    stopifnot(abs(sum(a$center * solve(v, a$center - d))) < 1e-10)
  }
}
# Independent linear-null oracle: endpoint roots must equal chi-square-1 cutoff.
for (d in list(c(4, 6), c(-4, 6), c(.4, 8))) {
  ci <- ws_interval(d, v)
  for (p in ci[ci > 0 & ci < 1]) {
    stats <- vapply(c(-1, 1), function(s) {
      a <- c(s * (1 - p)^2, -p^2)
      sum(a * d)^2 / drop(crossprod(a, v %*% a))
    }, numeric(1))
    stopifnot(abs(min(stats) - qchisq(.95, 1)) < 1e-7)
  }
  stopifnot(max(abs(ws_interval(d * 10, v * 100) - ci)) < 1e-9)
}
cat("Boundary selection, independent proxy precision and Fieller oracle checks passed.\n")
