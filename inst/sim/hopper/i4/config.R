# config.R -- simulation configuration for the P_med paper
# Mirrors the product-of-three config pattern: scenarios + run scales in one place.
# Sourced by run_simulation.R and analyze.R.

# ---- Run scales -------------------------------------------------------------
# test  : fast smoke check (seconds-minutes); NOT for reporting
# local : full grid on a workstation (slow; hours)
# cluster: full grid as a SLURM array (see submit_array.sh, to add)
SCALES <- list(
  test    = list(N = c(100, 250),                  reps = 50,   K = 2000,  n_boot = 200),
  local   = list(N = c(50, 100, 250, 500, 1000),   reps = 500,  K = 5000,  n_boot = 500),
  cluster = list(N = c(50, 100, 250, 500, 1000),   reps = 2000, K = 5000,  n_boot = 1000)
)

# Rep-chunks per (scenario, N) cell for the Hopper SLURM array (code/submit_array.sh):
# 7 scenarios x 5 N = 35 cells; 35 * N_CHUNKS = the array size.
N_CHUNKS <- 10

# ---- Data-generating mechanisms (DGMs) --------------------------------------
# Each scenario provides:
#   gen(N)    -> observed data.frame(X, M, Y)
#   po(N)     -> cross-world potential outcomes data.frame(Y1, Y0) under the TRUE
#                DGM, with Y1 = Y(1, M(1)) and Y0 = Y(1, M(0)) drawn independently
#                (the g-computation contrast). true_pmed() uses po(), so the truth
#                is the DGM estimand itself -- exact even where the estimation
#                model is misspecified (s5 mixture; logistic fit on probit binary).
#   family_y  -> "gaussian" or "binomial" (outcome working model for estimation)
#   label     -> human-readable name
SCENARIOS <- list(
  s1 = list(
    label = "Homogeneous (continuous Y)", family_y = "gaussian",
    gen = function(N) { X <- rbinom(N, 1, .5); M <- 0.5 * X + rnorm(N)
      Y <- 0.3 * X + 0.4 * M + rnorm(N); data.frame(X, M, Y) },
    po  = function(N) { M1 <- 0.5 + rnorm(N); M0 <- rnorm(N)
      data.frame(Y1 = 0.3 + 0.4 * M1 + rnorm(N), Y0 = 0.3 + 0.4 * M0 + rnorm(N)) }),
  s2 = list(
    label = "Inconsistent mediation (continuous Y)", family_y = "gaussian",
    gen = function(N) { X <- rbinom(N, 1, .5); M <- 0.5 * X + rnorm(N)
      Y <- -0.3 * X + 0.4 * M + rnorm(N); data.frame(X, M, Y) },
    po  = function(N) { M1 <- 0.5 + rnorm(N); M0 <- rnorm(N)
      data.frame(Y1 = -0.3 + 0.4 * M1 + rnorm(N), Y0 = -0.3 + 0.4 * M0 + rnorm(N)) }),
  s3 = list(
    # Intercept -2.06875 solved via uniroot() against the mixture-normal prevalence
    # formula (adversarial-audit finding #18, 2026-07-10): the prior -2.5 produced
    # ~2.1% prevalence, not the ~5% the prose claimed. Verified at N=2e6: 5.01%.
    label = "Binary Y, rare (~5%)", family_y = "binomial",
    gen = function(N) { X <- rbinom(N, 1, .5); M <- 0.5 * X + rnorm(N)
      ystar <- -2.06875 + 0.3 * X + 0.4 * M + rnorm(N)
      data.frame(X, M, Y = as.integer(ystar > 0)) },
    po  = function(N) { M1 <- 0.5 + rnorm(N); M0 <- rnorm(N)
      data.frame(Y1 = as.integer(-2.06875 + 0.3 + 0.4 * M1 + rnorm(N) > 0),
                 Y0 = as.integer(-2.06875 + 0.3 + 0.4 * M0 + rnorm(N) > 0)) }),
  s4 = list(
    # Intercept -0.25 is exact (not approximate): by symmetry of Phi, averaging
    # Phi((c+0.5X)/sd) over X in {0,1} equals 0.5 exactly when c = -0.25 (adversarial-
    # audit finding #18, 2026-07-10). Prior -0.5 produced ~41.1%, not ~50%.
    label = "Binary Y, balanced (~50%)", family_y = "binomial",
    gen = function(N) { X <- rbinom(N, 1, .5); M <- 0.5 * X + rnorm(N)
      ystar <- -0.25 + 0.3 * X + 0.4 * M + rnorm(N)
      data.frame(X, M, Y = as.integer(ystar > 0)) },
    po  = function(N) { M1 <- 0.5 + rnorm(N); M0 <- rnorm(N)
      data.frame(Y1 = as.integer(-0.25 + 0.3 + 0.4 * M1 + rnorm(N) > 0),
                 Y0 = as.integer(-0.25 + 0.3 + 0.4 * M0 + rnorm(N) > 0)) }),
  s5 = list(
    label = "Heterogeneous effects (mixture)", family_y = "gaussian",
    gen = function(N) { X <- rbinom(N, 1, .5); M <- 0.5 * X + rnorm(N)
      b <- ifelse(rbinom(N, 1, 0.7) == 1, 0.8, 0.0)
      Y <- 0.1 * X + b * M + rnorm(N); data.frame(X, M, Y) },
    po  = function(N) { M1 <- 0.5 + rnorm(N); M0 <- rnorm(N)
      # subgroup membership (responder slope) is a PERSON-level trait, so the
      # same b governs both cross-world arms Y(1,M(1)) and Y(1,M(0)).
      b  <- ifelse(rbinom(N, 1, 0.7) == 1, 0.8, 0.0)
      data.frame(Y1 = 0.1 + b * M1 + rnorm(N), Y0 = 0.1 + b * M0 + rnorm(N)) }),
  s6 = list(
    label = "Null mediation (b=0; no M->Y)", family_y = "gaussian",
    gen = function(N) { X <- rbinom(N, 1, .5); M <- 0.5 * X + rnorm(N)
      Y <- 0.3 * X + 0.0 * M + rnorm(N); data.frame(X, M, Y) },
    po  = function(N) { M1 <- 0.5 + rnorm(N); M0 <- rnorm(N)
      data.frame(Y1 = 0.3 + 0.0 * M1 + rnorm(N), Y0 = 0.3 + 0.0 * M0 + rnorm(N)) }),
  # s7: the OTHER composite-null branch, {a = 0}. X does not affect M (a = 0) but
  # M affects Y (b = 0.4). NIE = 0 and true P_med = 0.5 (M1, M0 i.i.d. -> the cross-world
  # difference D is symmetric about 0). This is the cell that exercises Branch 2 of the
  # MBCO two-branch constrained fit (mbco_pmed.R); pre-fix it had 0% null coverage.
  # NOT a manuscript Table-3 / Figure-2 scenario -- a diagnostic to confirm the a=0 fix at scale.
  s7 = list(
    label = "Null mediation (a=0; no X->M)", family_y = "gaussian",
    gen = function(N) { X <- rbinom(N, 1, .5); M <- 0.0 * X + rnorm(N)
      Y <- 0.3 * X + 0.4 * M + rnorm(N); data.frame(X, M, Y) },
    po  = function(N) { M1 <- rnorm(N); M0 <- rnorm(N)
      data.frame(Y1 = 0.3 + 0.4 * M1 + rnorm(N), Y0 = 0.3 + 0.4 * M0 + rnorm(N)) }),
  # s8: manuscript "Mechanism 7" (near-boundary), added post-adversarial-audit finding #23
  # (Phase 2 of docs/PLAN-adversarial-audit-fixes-2026-07-10.md, grill Branch 1). A large,
  # well-separated (a, b) combination in the same Gaussian single-mediator family as s1/s2/s5,
  # designed via the closed form (Appendix B) so true P_med lands in [0.85, 0.95] by
  # construction, not by guessing a "large effect": d = ab / sqrt(2(b^2 V_M + V_Y)) with
  # a = 2.5, b = 1.5, V_M = 1, V_Y = 2 gives Pmed = Phi(d) = 0.9008 (Wolfram-verified). Tests
  # whether MBCO's interior-coverage/width advantage over the percentile bootstrap holds near
  # the P_med boundary, the claim finding #23 flagged as unsupported by the six Table-3 mechanisms.
  s8 = list(
    label = "Near-boundary mediation (true Pmed ~ 0.90)", family_y = "gaussian",
    gen = function(N) { X <- rbinom(N, 1, .5); M <- 2.5 * X + rnorm(N)
      Y <- 0.3 * X + 1.5 * M + rnorm(N, sd = sqrt(2)); data.frame(X, M, Y) },
    po  = function(N) { M1 <- 2.5 + rnorm(N); M0 <- rnorm(N)
      data.frame(Y1 = 0.3 + 1.5 * M1 + rnorm(N, sd = sqrt(2)),
                 Y0 = 0.3 + 1.5 * M0 + rnorm(N, sd = sqrt(2))) })
)

OUT_DIR <- Sys.getenv("PMED_OUT_DIR", "simulation_results")   # truth.csv + all_results.csv written here
# Bootstrap Monte Carlo control: "crn" (common random numbers + control-variate
# centering, 2026-10-14 rerun) or "fresh" (original runs). See simulation_functions.R.
MC_CONTROL <- Sys.getenv("PMED_MC_CONTROL", "crn")
# Cross-world pairs per bootstrap-resample evaluation under "crn" (point estimate: 100 * K).
K_BOOT <- as.integer(Sys.getenv("PMED_K_BOOT", "50000"))
