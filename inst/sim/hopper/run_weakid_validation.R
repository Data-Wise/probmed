## #11 weak-ID / regularity THRESHOLD-VALIDATION grid (hopper).
##
## Purpose: move ward_residual()'s `weak_id_ratio_threshold` (default 3) and
## `oe_snr_threshold` (default 2) from PROVISIONAL to validated (or corrected).
##
## OUTCOME (ran 2026-07-16, job 4277259, 48k reps; "corrected" is what happened):
## A2's coverage separation is real in all 24 cells but ~97% width-tautological;
## A1 flags OVER-covering (uninformatively wide) intervals and is the stronger
## gate; t in 1.75-4 indistinguishable, 3 kept as convention. See
## ../README.md "What the grid found" and ../results/weakid_validation_cells.csv.
##
## Design notes (why this grid differs from run_gauge_boot_grid.R):
##  * Axis. The gaugeboot grid varies n x tint x binY, which lands oe_snr at the
##    EXTREMES (~1 or ~12). A threshold arbitrates in the MIDDLE (oe_snr ~1.5-5),
##    which no gaugeboot cell samples. So this grid sweeps the A-effect SCALE `s`,
##    the axis that drives oe_snr monotonically through the crossover.
##  * Recorded quantities. We record the SHIPPED fields (weak_id, weak_id_ratio,
##    oe_snr, oe_regular) straight off the fitted object, so we validate exactly
##    what users get -- not a re-derived proxy. (The gaugeboot grid predates these
##    fields and recorded an SE ratio + max-endpoint divergence instead, which is
##    why it could only corroborate the threshold indirectly.)
##  * B = 200 = ward_residual()'s DEFAULT. A user-facing flag must be calibrated
##    under the conditions users actually run. (Probe: B=999 raises the ratio only
##    ~6-8% vs B=200, so the threshold is not strongly B-dependent.)
##
## Analysis target (see collate_weakid_validation.R): percentile coverage is ~1.00
## everywhere, so the flag CANNOT be validated against percentile under-coverage.
## The flag's actual claim is "W is weakly identified; its CI is least trustworthy
## here", and the CI users get BY DEFAULT is the analytic Wald one (~0.86-0.91).
## So the validation is an operating-characteristic question:
##     Is Wald coverage materially worse when weak_id = TRUE than when FALSE?
## A validated threshold => unflagged group ~ nominal, flagged group clearly worse.
suppressMessages(library(probmed))

## --- stale-package guard: the shipped weak-ID fields must EXIST and bootstrap must be live.
## Without this a wrong R_LIBS silently yields NA flags (old package) or a no-op bootstrap,
## and the whole grid would be measuring nothing. Fail loudly, never silently ship.
.chk <- local({
  set.seed(99); n <- 300; C <- rnorm(n); A <- rbinom(n, 1, plogis(0.3 * C))
  M <- 0.5 * A + 0.3 * C + rnorm(n); Y <- 0.4 * A + 0.6 * M + rnorm(n)
  d <- data.frame(A, M, Y, C)
  b <- tryCatch(ward_residual(d, se_method = "bootstrap", B = 50L), error = function(e) NULL)
  a <- tryCatch(ward_residual(d, se_method = "analytic"), error = function(e) NULL)
  if (is.null(a) || is.null(b)) return(FALSE)
  has_fields <- all(c("weak_id", "weak_id_ratio", "oe_snr", "oe_regular") %in%
                      S7::prop_names(b))
  live_boot <- !isTRUE(all.equal(a@W_ci, b@W_ci))
  has_fields && live_boot && !is.na(b@weak_id_ratio)
})
if (!isTRUE(.chk))
  stop("STALE probmed: weak-ID fields missing/NA or bootstrap is a no-op. ",
       "Needs dev >= PR #23 installed FIRST on R_LIBS.")

expit <- function(x) 1 / (1 + exp(-x))

## DGP: `s` scales BOTH A->M and A->Y (incl. the A*M interaction), sweeping
## identification strength from near-null (s~0) to strong (s=1).
gen <- function(n, s, binY) {
  C <- rnorm(n); A <- rbinom(n, 1, expit(-0.2 + 0.8 * C))
  M <- s * 0.6 * A + 0.4 * C + rnorm(n)
  lin <- s * 0.5 * A + 0.7 * M + s * 0.4 * A * M + 0.3 * C
  Y <- if (binY) rbinom(n, 1, expit(lin)) else lin + rnorm(n)
  data.frame(A, M, Y, C)
}

## EXACT truth -- no Monte Carlo. `lin` is a linear combination of normals, so it
## is exactly Gaussian even in the binary case:
##   lin = 0.5s*a + b*M + 0.3C,   b = 0.7 + 0.4*s*a,   M = 0.6s*ap + 0.4C + eM
##       = 0.5s*a + 0.6s*b*ap + (0.4b + 0.3)C + b*eM,  C, eM ~ iid N(0,1)
##   => lin ~ N(mu, sd^2),  mu = 0.5s*a + 0.6s*b*ap,  var = (0.4b+0.3)^2 + b^2
## Continuous Y: theta = mu. Binary Y: theta = E[expit(lin)] -- a 1-D Gaussian
## integral, by quadrature to ~1e-10.
##
## WHY NOT MONTE CARLO (do not reintroduce it): an earlier version drew N = 2e6
## per corner. W = R/OE divides by OE -> 0, so MC error in the denominator
## explodes the ratio exactly where this grid is trying to measure -- and the
## call was unseeded, so each of the 8 chunks in a cell scored coverage against
## a DIFFERENT draw of the truth (the collator then reported only the first).
## Scale of the problem at s = 0.05, measured over repeated MC runs: the MC
## truth has SD ~0.017-0.019 and ranges roughly [-0.015, +0.045] -- larger than
## the estimand itself (exact: 0.0129 continuous / 0.0101 binary) and unstable
## in SIGN. (One MC run returned 0.043 vs exact 0.0129; quoting that single draw
## as "the" MC value would repeat exactly the single-draw error this comment
## warns about -- the point is the spread, not any one number.) MC converges
## only at large s, e.g. s = 0.5: MC ~0.118 vs exact 0.1154.
##
## NOTE the sibling grid (run_gauge_boot_grid.R) legitimately keeps MC: its axis
## is `tint`, so OE = 0.92 + 0.6*tint is never near 0, and its stored trW is off
## by only ~0.002 (<= 0.04 of empSD_W). The hazard is specific to sweeping
## toward the null.
theta_exact <- function(a, ap, s, binY) {
  b  <- 0.7 + 0.4 * s * a
  mu <- 0.5 * s * a + 0.6 * s * b * ap
  sd <- sqrt((0.4 * b + 0.3)^2 + b^2)
  if (!binY) return(mu)
  stats::integrate(function(x) expit(x) * stats::dnorm(x, mu, sd),
                   lower = mu - 12 * sd, upper = mu + 12 * sd,
                   rel.tol = 1e-10)$value
}
truth <- function(s, binY) {
  t11 <- theta_exact(1, 1, s, binY); t10 <- theta_exact(1, 0, s, binY)
  t01 <- theta_exact(0, 1, s, binY); t00 <- theta_exact(0, 0, s, binY)
  OE <- t11 - t00
  c(W = (OE - (t10 - t00) - (t01 - t00)) / OE, P = (t01 - t00) / OE)
}

## 24 cells x 8 chunks = 192 array tasks; 250 reps/chunk => nsim = 2000 per cell.
cells <- expand.grid(s = c(0.05, 0.10, 0.15, 0.20, 0.30, 0.50),
                     n = c(800, 3000), binY = c(FALSE, TRUE))
NCHUNK <- 8L; REPS_PER <- 250L; B <- 200L    # B = shipped default

aid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))   # 1..192
ci  <- ((aid - 1L) %/% NCHUNK) + 1L                          # cell  1..24
ch  <- ((aid - 1L) %%  NCHUNK) + 1L                          # chunk 1..8
s <- cells$s[ci]; n <- cells$n[ci]; binY <- cells$binY[ci]
tr <- truth(s, binY)

base <- 100000L * ci + 1000L * ch
rows <- vector("list", REPS_PER); nfail <- 0L
for (j in seq_len(REPS_PER)) {
  seed <- base + j
  set.seed(seed); d <- gen(n, s, binY)
  fa <- tryCatch(ward_residual(d, seed = seed), error = function(e) NULL)
  fb <- tryCatch(suppressWarnings(
          ward_residual(d, seed = seed, se_method = "bootstrap", B = B)),
        error = function(e) NULL)
  if (is.null(fa) || is.null(fb)) { nfail <- nfail + 1L; next }
  rows[[j]] <- data.frame(
    cell = ci, chunk = ch, seed = seed, n = n, s = s, binY = binY,
    trW = tr["W"], W = fa@W,
    ## --- the SHIPPED diagnostic fields, verbatim ---
    weak_id       = fb@weak_id,
    weak_id_ratio = fb@weak_id_ratio,
    oe_snr        = fb@oe_snr,
    oe_regular    = fb@oe_regular,
    ## --- coverage of the DEFAULT (Wald) interval and the percentile one ---
    covW_an  = tr["W"] >= fa@W_ci[1] && tr["W"] <= fa@W_ci[2],
    covW_pct = tr["W"] >= fb@W_ci[1] && tr["W"] <= fb@W_ci[2],
    ## widths, so a threshold sweep can be redone without re-running the grid
    wid_wald = diff(fb@W_ci_wald), wid_pct = diff(fb@W_ci),
    row.names = NULL)
}
out <- do.call(rbind, rows)
attr(out, "nfail") <- nfail
dir.create("~/weakid_val/parts", recursive = TRUE, showWarnings = FALSE)
f <- sprintf("~/weakid_val/parts/part_cell%02d_chunk%02d.rds", ci, ch)
saveRDS(out, f)
cat(sprintf("cell %d (s=%.2f n=%d binY=%s) chunk %d: %d rows, nfail=%d -> %s\n",
            ci, s, n, binY, ch, nrow(out), nfail, f))
