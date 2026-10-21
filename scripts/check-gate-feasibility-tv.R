#!/usr/bin/env Rscript
# T3 feasibility check: can ANY [0,1]-valued estimator pass the frozen bias/SD < .10 gate in
# both cell 5 (direct = 0, truth P = 1) and cell 7 (direct = .05, truth P = .7721) at n = 500?
#
# POPULATION argument (two-point / total variation), for a possibly randomized estimator T in [0,1]:
#   cell 5: |b5| < .1 SD5 and SD5 <= sqrt(E5 (1 - E5)), b5 = 1 - E5   =>  E5 >= .9901.
#   cell 7: |E7 - .7721| < .1 SD7, SD7 <= sqrt(E7 (1 - E7))           =>  E7 <= .811 (loose: .8221).
#   |E5 - E7| <= TV(P5, P7)  =>  passing both needs TV >= .9901 - .811 = .179 (loose: .168).
# EMPIRICAL gate (R = 2000 replicates, sample SD with n-1): for values in [0,1] the sample SD
#   satisfies s^2 <= R/(R-1) * m (1 - m) with m the sample mean, so the same algebra holds for
#   the sample means; passing both empirical gates needs mbar5 - mbar7 >= .179 (approximately,
#   the R/(R-1) factor is 1.0005). The population gap is E5 - E7 <= TV, so an empirical pass needs
#   sampling noise in (mbar5 - mbar7) of at least .179 - TV. Bernstein bound with worst-case
#   variance is reported below.
# The two data laws differ only in the Y | A, M, C logit (shift .05 * A), so
#   TV = E_5[(1 - Q/P)_+], log(Q/P) = sum of Bernoulli log-likelihood ratios (MC, arbitrary seed 97).
# Controls: shift 0 must give TV = 0; shifts .25 and .5 must give TV well above .179.
set.seed(97)
n <- 500L; R <- 20000L; alpha <- .7
tv_shift <- function(shift) {
  v <- vapply(seq_len(R), function(i) {
    C <- stats::rnorm(n); A <- stats::rbinom(n, 1, stats::plogis(.3 * C))
    M <- alpha * A + .5 * C + stats::rnorm(n)
    e5 <- -.5 + .8 * M + .3 * C
    e7 <- -.5 + shift * A + .8 * M + .3 * C
    Y <- stats::rbinom(n, 1, stats::plogis(e5))
    ll <- sum(stats::dbinom(Y, 1, stats::plogis(e7), log = TRUE) -
      stats::dbinom(Y, 1, stats::plogis(e5), log = TRUE))
    max(0, 1 - exp(ll))
  }, 0)
  c(tv = mean(v), se = sd(v) / sqrt(R))
}
res <- t(vapply(c(0, .05, .25, .5), tv_shift, c(tv = 0, se = 0)))
res <- data.frame(shift = c(0, .05, .25, .5), res)
cat("TV(cell 5 data law, shifted data law) at n = 500, R = 20000 datasets each:\n")
print(res, row.names = FALSE, digits = 4)
tv <- res$tv[res$shift == .05]; se <- res$se[res$shift == .05]
cat(sprintf("\nControl shift 0 gives TV = %.4f (must be 0): %s\n", res$tv[1], if (res$tv[1] == 0) "ok" else "BAD"))
cat(sprintf("Controls .25 and .5 exceed the needed .179: %s\n",
  if (all(res$tv[res$shift %in% c(.25, .5)] > .179)) "ok" else "BAD"))

need_pop <- .9901 - .811
cat(sprintf("\nPOPULATION: TV(cell 5, cell 7) = %.4f (MC SE %.4f); needed >= %.4f (loose %.4f). %s\n",
  tv, se, need_pop, .9901 - .8221,
  if (tv + 4 * se < need_pop) "Population gate infeasible for any [0,1] estimator." else "Not ruled out."))

# EMPIRICAL: an empirical pass needs BOTH sample means to satisfy their one-sided conditions
# separately: mbar5 >= .9901 and mbar7 <= .8112 (cell 7 also needs a lower condition, ignored: this
# only enlarges the bound). The cells use independent datasets, so P(pass) <= P5 * P7 with
# P5 = P(mbar5 >= c5), P7 = P(mbar7 <= c7). The population means must satisfy |E5 - E7| <= TV.
# Worst case over E5 (E7 = E5 - TV, the value most favorable to a pass), using Bernstein
# (variance E(1-E)/R, range 1/R) and the normal approximation. Replaces an earlier difference-of-
# means argument whose worst-case variance was mis-specified (independent audit, 2026-10-21).
Rr <- 2000; c5 <- .9901; c7 <- .8112; tvw <- tv
bern1 <- function(t, v) if (t <= 0) 1 else exp(-t^2 / (2 * (v + (1 / Rr) * t / 3)))
norm1 <- function(t, v) if (t <= 0) 1 else stats::pnorm(-t / sqrt(v))
grid <- seq(c7 + .0005, 1, by = .0005)
best <- c(bern = 0, norm = 0, e5_bern = NA, e5_norm = NA)
for (e5 in grid) {
  e7 <- e5 - tvw
  if (e7 <= 0 || e7 >= 1) next
  v5 <- e5 * (1 - e5) / Rr; v7 <- e7 * (1 - e7) / Rr
  pb <- bern1(c5 - e5, v5) * bern1(e7 - c7, v7)
  pn <- norm1(c5 - e5, v5) * norm1(e7 - c7, v7)
  if (pb > best["bern"]) { best["bern"] <- pb; best["e5_bern"] <- e5 }
  if (pn > best["norm"]) { best["norm"] <- pn; best["e5_norm"] <- e5 }
}
cat(sprintf("EMPIRICAL (R = %d independent replicates per cell): P(empirical pass) <= P(mbar5 >= %.4f) * P(mbar7 <= %.4f)\n",
  Rr, c5, c7))
cat(sprintf("  worst case over population means (E7 = E5 - TV): Bernstein <= %.2e (at E5 = %.4f); normal approximation <= %.2e (at E5 = %.4f)\n",
  best["bern"], best["e5_bern"], best["norm"], best["e5_norm"]))
cat("  Applies to any single estimator with values in [0,1] evaluated on independent datasets in the two cells.\n")
