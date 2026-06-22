## Tests for incr_pmed() / IncrPmedResult (incremental mediated elasticity, Paper 2)
.ip_gen <- function(n, tint, seed = 1) {
  set.seed(seed)
  C <- rnorm(n)
  A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
  M <- 0.6 * A + 0.4 * C + rnorm(n)
  Y <- 0.5 * A + 0.7 * M + tint * A * M + 0.3 * C + rnorm(n)
  data.frame(A, M, Y, C)
}

# Plug-in oracle for the headline properties: nested corner means nu(a, a', C),
# the elasticity decomposition, and Theta(d, d) for a numeric derivative check.
# This mirrors the validated prototype (sims/incr_pmed_proto.R) and lets us test
# the exact-additivity (Theorem 1) and reduction (Theorem 2) theorems directly,
# rather than asserting the tautology tot == dir + med on the one-step output.
.ip_oracle <- function(d, deltas) {
  g <- stats::predict(stats::glm(A ~ C, data = d, family = stats::binomial()),
                      type = "response")
  om <- stats::glm(Y ~ A * M + C, data = d, family = stats::gaussian())
  nu <- matrix(NA_real_, nrow(d), 4, dimnames = list(NULL, c("11", "10", "01", "00")))
  for (cc in list(c(1, 1, "11"), c(1, 0, "10"), c(0, 1, "01"), c(0, 0, "00"))) {
    a <- as.numeric(cc[1]); ap <- as.numeric(cc[2]); nm <- cc[3]
    muAM <- stats::predict(om, newdata = transform(d, A = a))
    etam <- stats::lm(yy ~ C, data = data.frame(yy = muAM[d$A == ap], C = d$C[d$A == ap]))
    nu[, nm] <- stats::predict(etam, newdata = d)
  }
  Theta <- function(del) {
    q <- del * g / (del * g + 1 - g)
    mean(q * q * nu[, "11"] + q * (1 - q) * nu[, "10"] +
         (1 - q) * q * nu[, "01"] + (1 - q) * (1 - q) * nu[, "00"])
  }
  t(vapply(deltas, function(del) {
    q  <- del * g / (del * g + 1 - g)
    qp <- g * (1 - g) / (del * g + 1 - g)^2
    dir <- mean(qp * (q * (nu[, "11"] - nu[, "01"]) + (1 - q) * (nu[, "10"] - nu[, "00"])))
    med <- mean(qp * (q * (nu[, "11"] - nu[, "10"]) + (1 - q) * (nu[, "01"] - nu[, "00"])))
    tot_num <- (Theta(del + 1e-5) - Theta(del - 1e-5)) / 2e-5
    c(delta = del, dir = dir, med = med, tot_analytic = dir + med,
      tot_numeric = tot_num, remainder = tot_num - (dir + med), Pmed = med / (dir + med))
  }, numeric(7)))
}

test_that("incr_pmed returns an IncrPmedResult with a finite curve", {
  r <- incr_pmed(.ip_gen(1500, 0.0), deltas = c(0.5, 1, 2))
  expect_true(S7::S7_inherits(r, IncrPmedResult))
  expect_equal(nrow(r@curve), 3)
  expect_true(all(is.finite(r@curve$Pmed)))
  expect_true(all(is.finite(r@curve$se)))
  expect_true(all(is.finite(r@curve$lo)) && all(is.finite(r@curve$hi)))
  expect_true(all(r@curve$lo <= r@curve$hi))
})

test_that("Theorem 1: exact additivity, remainder R(delta) = 0 for every delta", {
  # numeric d/ddelta Theta(delta, delta) equals the analytic dir + med
  o <- .ip_oracle(.ip_gen(20000, 0.9), deltas = c(1 / 3, 1 / 2, 1, 2, 3))
  expect_true(all(abs(o[, "remainder"]) < 1e-3))
  expect_equal(o[, "tot_numeric"], o[, "tot_analytic"], tolerance = 1e-3,
               ignore_attr = TRUE)
})

test_that("Theorem 2: reduction to classical P_med under no A*M interaction", {
  # linear-Gaussian, no interaction => flat curve at b_a * t_m / (t_a + b_a * t_m)
  classical <- 0.6 * 0.7 / (0.5 + 0.6 * 0.7)   # = 0.4565
  o <- .ip_oracle(.ip_gen(40000, 0.0), deltas = c(1 / 3, 1 / 2, 1, 2, 3))
  expect_true(all(abs(o[, "Pmed"] - classical) < 0.03))  # at the classical value
  expect_lt(diff(range(o[, "Pmed"])), 0.01)              # and flat across delta
})

test_that("curve bends with delta when an A*M interaction is present", {
  o <- .ip_oracle(.ip_gen(40000, 0.9), deltas = c(1 / 3, 1, 3))
  expect_gt(diff(range(o[, "Pmed"])), 0.02)
})

test_that("one-step estimator recovers the classical value at delta = 1 (no interaction)", {
  classical <- 0.6 * 0.7 / (0.5 + 0.6 * 0.7)
  r <- incr_pmed(.ip_gen(8000, 0.0), deltas = 1, K = 5L)
  expect_true(abs(r@curve$Pmed[1] - classical) < 0.08)
  expect_lte(r@curve$lo[1], classical)
  expect_gte(r@curve$hi[1], classical)
})

# ---- Estimator-level coverage of the headline properties (remediation A-3) ----
# The Theorem 1/2 tests above run on .ip_oracle(); these exercise the SHIPPED
# incr_pmed() so the headline claims are verified on the estimator users call.

test_that("Theorem 1 (estimator): returned curve satisfies dir + med == tot", {
  # NOTE: in incr_pmed() tot is *defined* as dir + med (R/incremental-pmed.R),
  # so on the shipped @curve this identity holds to machine precision -- it is a
  # structural consistency check, not an empirical test of the chain rule. The
  # genuine, non-tautological additivity test (numeric d/ddelta Theta == dir+med)
  # remains the oracle test above, which is the only way to validate Theorem 1
  # without re-exposing Theta(delta, delta) from the estimator internals.
  r <- incr_pmed(.ip_gen(12000, 0.9), deltas = c(0.5, 1, 2), K = 5L)
  expect_equal(r@curve$dir + r@curve$med, r@curve$tot, tolerance = 1e-12)
  expect_equal(r@curve$Pmed, r@curve$med / r@curve$tot, tolerance = 1e-12)
})

test_that("Theorem 2 (estimator): no-interaction curve is flat at classical P_med", {
  # large-n, no A*M interaction => incr_pmed()'s own curve should be ~flat and
  # ~ the classical proportion mediated (not just the oracle's).
  classical <- 0.6 * 0.7 / (0.5 + 0.6 * 0.7)        # = 0.4565
  r <- incr_pmed(.ip_gen(18000, 0.0), deltas = c(0.5, 1, 2), K = 5L)
  expect_lt(max(abs(r@curve$Pmed - classical)), 0.06)  # near classical value
  expect_lt(diff(range(r@curve$Pmed)), 0.05)           # and flat across delta
})

test_that("curve bends (estimator): strong interaction => curve range exceeds threshold", {
  # under a strong A*M interaction the estimator's own curve must vary with delta;
  # 0.05 is the print() "bends with delta" threshold.
  r <- incr_pmed(.ip_gen(18000, 0.9), deltas = c(1 / 3, 1, 3), K = 5L)
  expect_gt(diff(range(r@curve$Pmed)), 0.05)
})

test_that("multiple covariates are supported", {
  d <- .ip_gen(1500, 0.3); d$C2 <- rnorm(nrow(d))
  r <- incr_pmed(d, deltas = c(0.5, 2), covars = c("C", "C2"))
  expect_true(all(is.finite(r@curve$Pmed)))
})

test_that("print method runs and returns invisibly", {
  r <- incr_pmed(.ip_gen(1000, 0.0), deltas = c(0.5, 1, 2))
  expect_output(print(r), "Incremental mediated elasticity")
  expect_identical(withVisible(print(r))$visible, FALSE)
})

# ---- Kennedy (2019) Cor. 2 g-score term validation ----

test_that("g-score correction: internal terms are mean-zero (Neyman orthogonality)", {
  # E[(dq/dg)*(A - g(C))*gamma(C)] = 0 when g is correctly specified.
  # With correctly specified g and large n, the t-stat should be < 3.5.
  skip_on_cran()
  set.seed(99); n <- 8000L
  d <- .ip_gen(n, 0.0)
  fit <- probmed:::.corner_fit(d, K = 5L, binY = FALSE, covars = "C")
  phi <- fit$phi; g <- fit$g; del <- 1.0
  q    <- del * g / (del * g + 1 - g)
  dqg  <- del / (del * g + 1 - g)^2
  resid <- d$A - g
  gamma_dir <- q * (phi[, "11"] - phi[, "01"]) + (1 - q) * (phi[, "10"] - phi[, "00"])
  gamma_med <- q * (phi[, "11"] - phi[, "10"]) + (1 - q) * (phi[, "01"] - phi[, "00"])
  gterm_dir <- dqg * gamma_dir * resid
  gterm_med <- dqg * gamma_med * resid
  t_dir <- abs(mean(gterm_dir)) * sqrt(n) / sd(gterm_dir)
  t_med <- abs(mean(gterm_med)) * sqrt(n) / sd(gterm_med)
  expect_lt(t_dir, 3.5)
  expect_lt(t_med, 3.5)
})

test_that("g-score EIF: CI covers oracle truth at delta=1 in large-n no-interaction DGP", {
  skip_on_cran()
  classical <- 0.6 * 0.7 / (0.5 + 0.6 * 0.7)
  r <- incr_pmed(.ip_gen(12000, 0.0), deltas = 1, K = 5L)
  expect_gte(r@curve$hi[1], classical)
  expect_lte(r@curve$lo[1], classical)
})

test_that("g-score EIF: SE is positive and plausible after adding g-score correction", {
  skip_on_cran()
  r <- incr_pmed(.ip_gen(6000, 0.4), deltas = c(0.5, 1, 2), K = 5L)
  expect_true(all(is.finite(r@curve$se) & r@curve$se > 0))
  expect_true(all(r@curve$se < 0.3))  # implausibly large SE would signal a bug
})

test_that("g-score EIF: nominal 95% CI achieves coverage >= 80% in small simulation", {
  skip_on_cran()
  classical <- 0.6 * 0.7 / (0.5 + 0.6 * 0.7)
  set.seed(77)
  seeds <- sample.int(1e5, 25)
  covered <- vapply(seeds, function(s) {
    r <- incr_pmed(.ip_gen(2500, 0.0, seed = s), deltas = 1, K = 5L)
    r@curve$lo[1] <= classical && classical <= r@curve$hi[1]
  }, logical(1))
  expect_gte(mean(covered), 0.80)
})
