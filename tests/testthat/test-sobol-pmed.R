## Tests for sobol_pmed() / SobolPmedResult (Sobol / variance-scale P_med^{sigma2}).
## Adapted from the verified sim suite (03-sobol-pmed/sims/test-sobol-pmed.R); these
## are SELF-CONTAINED (no external source(): DGP, cells and closed-form truth are
## inlined, mirroring test-gauge-pmed.R) and call the PACKAGE generic sobol_pmed().

.sp_expit <- function(x) 1 / (1 + exp(-x))

## linear-Gaussian-with-interaction DGP; corner means (and thus the Sobol truth)
## are closed form.
.sp_gen <- function(n, p) with(p, {
  C <- rnorm(n); A <- rbinom(n, 1, .sp_expit(a0 + a1 * C))
  M <- b0 + ba * A + bc * C + rnorm(n, 0, sM)
  Y <- t0 + ta * A + tm * M + kappa * A * M + tc * C + rnorm(n, 0, sY)
  data.frame(C, A, M, Y)
})

## closed-form corner means theta(a, a') and the Sobol truth via the package helper.
.sp_theta <- function(p) with(p,
  c("11" = t0 + ta * 1 + (tm + kappa * 1) * (b0 + ba * 1),
    "10" = t0 + ta * 1 + (tm + kappa * 1) * (b0 + ba * 0),
    "01" = t0 + ta * 0 + (tm + kappa * 0) * (b0 + ba * 1),
    "00" = t0 + ta * 0 + (tm + kappa * 0) * (b0 + ba * 0)))
.sp_truth <- function(p, pd = .5, pm = .5) unname(sobol_from_theta(.sp_theta(p), pd, pm)["Pmed_sobol"])

.sp_base <- list(a0 = -0.2, a1 = 0.8, b0 = 0, bc = 0.4, t0 = 0, tc = 0.3, sM = 1, sY = 1)
cell_additive    <- modifyList(.sp_base, list(ba = .6, ta = .5, tm = .7, kappa = 0))   # no interaction
cell_interaction <- modifyList(.sp_base, list(ba = .6, ta = .5, tm = .4, kappa = .6))  # strong med + int
cell_null        <- modifyList(.sp_base, list(ba = .6, ta = .5, tm = 0,  kappa = 0))   # Dm=0 -> V_med=0

## -----------------------------------------------------------------------------------------------
test_that("returns a SobolPmedResult; P_med in [0,1] and variance components >= 0", {
  for (p in list(cell_additive, cell_interaction, cell_null)) {
    set.seed(11)
    fit <- sobol_pmed(.sp_gen(2000, p), seed = 1L, warn_boundary = FALSE)
    expect_true(S7::S7_inherits(fit, SobolPmedResult))
    expect_gte(fit@p_med, -1e-8); expect_lte(fit@p_med, 1 + 1e-8)
    expect_gte(fit@Vd,  -1e-8); expect_gte(fit@Vm, -1e-8); expect_gte(fit@Vdm, -1e-8)
    expect_gt(fit@VT, 0)
    expect_gte(fit@S1_med, -1e-8); expect_lte(fit@S1_med, 1 + 1e-8)
    expect_gte(fit@ST_med, -1e-8); expect_lte(fit@ST_med, 1 + 1e-8)
  }
})

## -----------------------------------------------------------------------------------------------
test_that("reduction (kappa=0): P_med_sobol ~ NIE^2 / (NIE^2 + NDE^2) at large n", {
  truth <- .sp_truth(cell_additive, pd = .5, pm = .5)
  NIE <- cell_additive$tm * cell_additive$ba   # 0.7 * 0.6
  NDE <- cell_additive$ta                       # 0.5
  reduction <- NIE^2 / (NIE^2 + NDE^2)
  ## the closed-form estimand truth must equal the analytic reduction at kappa=0, pd=pm=.5
  expect_equal(truth, reduction, tolerance = 1e-8)
  ## average seeds at n=4000 to beat single-seed MC error and DISCRIMINATE the Sobol
  ## truth (~0.414) from the classic share NIE/(NIE+NDE) (~0.457).
  ests <- vapply(1:10, function(s) {
    set.seed(100 + s)
    sobol_pmed(.sp_gen(4000, cell_additive), seed = s, warn_boundary = FALSE)@p_med
  }, numeric(1))
  expect_equal(mean(ests), truth, tolerance = 0.03)
})

## -----------------------------------------------------------------------------------------------
test_that("finite inference: se and ci finite, ci brackets the point estimate", {
  set.seed(21)
  fit <- sobol_pmed(.sp_gen(3000, cell_interaction), seed = 3L, warn_boundary = FALSE)
  expect_true(is.finite(fit@se)); expect_gt(fit@se, 0)
  expect_length(fit@ci, 2)
  expect_true(all(is.finite(fit@ci)))
  expect_lte(fit@ci[1], fit@p_med)
  expect_gte(fit@ci[2], fit@p_med)
})

## -----------------------------------------------------------------------------------------------
test_that("split-guard: rejects H0 V_med=0 under strong mediation; rarely under true null", {
  ## strong-mediation DGP (Dm ~ 0.42 >> 0): must reject H0: V_med = 0
  set.seed(31)
  fit_strong <- sobol_pmed(.sp_gen(2000, cell_interaction), seed = 5L, warn_boundary = FALSE)
  expect_identical(fit_strong@vmed_split_reject, 1L)
  ## true-null DGP (tm=kappa=0 => closed-form Dm=0 => V_med=0): split test is
  ## conservative (verified ~0.00); threshold 0.2 separates the fix from the old
  ## plug-in over-rejection (~0.25-0.43).
  rej <- vapply(1:20, function(s) {
    set.seed(300 + s)
    sobol_pmed(.sp_gen(1500, cell_null), seed = s, warn_boundary = FALSE)@vmed_split_reject
  }, integer(1))
  expect_lt(mean(rej), 0.2)
})

## -----------------------------------------------------------------------------------------------
test_that("boundary gating: Procedure-A upper bound reported when test does not reject", {
  set.seed(51)
  fit <- sobol_pmed(.sp_gen(1500, cell_null), seed = 1L, warn_boundary = FALSE,
                    procedure = "A")                 # A-specific test: report the gated interval
  if (isTRUE(fit@boundary)) {
    expect_identical(fit@ci[1], 0)
    expect_equal(fit@ci[2], fit@Pmed_upper)
    expect_gte(fit@Pmed_upper, 0)
  } else {
    ## off-boundary the reported CI equals the Wald CI
    expect_equal(fit@ci, fit@ci_wald)
  }
})

## -----------------------------------------------------------------------------------------------
test_that("warn_boundary=TRUE emits a warning at the boundary", {
  ## find a seed at the boundary, then confirm the warning fires there.
  set.seed(61)
  d <- .sp_gen(1500, cell_null)
  silent <- sobol_pmed(d, seed = 1L, warn_boundary = FALSE)
  if (isTRUE(silent@boundary)) {
    ## A path warns (non-uniform gated bound); B path (default) is uniform -> message, not warning
    expect_warning(sobol_pmed(d, seed = 1L, warn_boundary = TRUE, procedure = "A"), "Procedure A")
    expect_no_warning(sobol_pmed(d, seed = 1L, warn_boundary = TRUE, procedure = "B"))
  } else {
    succeed("seed not at boundary; warning path covered elsewhere")
  }
})

## -----------------------------------------------------------------------------------------------
test_that("Procedure B (default) returns the mapped split-CI image, bounded and gate-free", {
  set.seed(51)
  fit <- sobol_pmed(.sp_gen(1500, cell_null), seed = 1L, warn_boundary = FALSE)  # default procedure="B"
  expect_identical(fit@procedure, "B")
  expect_identical(fit@ci, fit@ci_B1)              # default reports the coherent same-sample image
  expect_true(all(is.finite(fit@ci)))
  expect_gte(fit@ci[1], 0)                          # image of a square is non-negative
  expect_lte(fit@ci[1], fit@ci[2])
  expect_gte(fit@ci[2], fit@p_med)                  # B1 always brackets the point estimate (centered at full Dm)
  ## B2 (decorrelated split) present as a robustness diagnostic
  expect_true(all(is.finite(fit@ci_B2)) && fit@ci_B2[1] >= 0 && fit@ci_B2[1] <= fit@ci_B2[2])
  ## B does NOT gate: its interval is the squared-map image, distinct from A's [0, Pmed_upper] rule
  expect_false(identical(fit@ci, fit@ci_A) && !isTRUE(fit@boundary))
})

## -----------------------------------------------------------------------------------------------
test_that("determinism: same seed -> identical point estimate and inference", {
  set.seed(41)
  d <- .sp_gen(2000, cell_interaction)           # generate ONCE; reuse for both calls
  f1 <- sobol_pmed(d, seed = 7L, warn_boundary = FALSE)
  f2 <- sobol_pmed(d, seed = 7L, warn_boundary = FALSE)
  expect_identical(f1@p_med, f2@p_med)
  expect_identical(f1@se, f2@se)
  expect_identical(f1@ci, f2@ci)
  expect_identical(f1@vmed_split_reject, f2@vmed_split_reject)
})

## -----------------------------------------------------------------------------------------------
test_that("sobol_from_theta matches the variance decomposition identities", {
  th <- .sp_theta(cell_interaction)
  v <- sobol_from_theta(th, pd = .5, pm = .5)
  expect_equal(unname(v["VT"]), unname(v["Vd"] + v["Vm"] + v["Vdm"]), tolerance = 1e-12)
  expect_equal(unname(v["Pmed_sobol"]), unname(v["Vm"] / v["VT"]), tolerance = 1e-12)
  expect_equal(unname(v["ST_med"]), unname((v["Vm"] + v["Vdm"]) / v["VT"]), tolerance = 1e-12)
})
