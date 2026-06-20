## wasserstein_pmed.R — transport-scale proportion mediated P_med^W (PROTOTYPE).
## pmed-modern Paper 4. Drop-in candidate for probmed/R/ (matches S7 conventions),
## NOT yet wired into NAMESPACE/DESCRIPTION — promote deliberately (probmed CRAN-blocked on medfit).
##
## Estimand (corner LAWS nu(a,a') = Law(Y(a, M(a'))), a=1, a*=0):
##   NIE^W = W2(nu_11, nu_10),  NDE^W = W2(nu_10, nu_00),  TE^W = W2(nu_11, nu_00)
##   P_med^W = NIE^W / (NIE^W + NDE^W) in [0,1];  synergy G^W = NIE^W + NDE^W - TE^W >= 0.
##
## STATUS (2026-06-19): transport + EIF engine VERIFIED (manuscript repo
##   04-wasserstein-pmed/sims/wasserstein_verify.R FACT1-6). d=1 only.
##   Corner-law estimator = Monte-Carlo g-computation (real, runnable). SE: nonparametric
##   bootstrap (idiomatic in probmed) + a verified closed-/numeric-form EIF SE for the
##   transport layer (known corner laws).
## TODO(Claude Code, R/ promotion): (a) cross-fit triply-robust corner-LAW EIF SE (reuse the
##   .gauge_fit pattern, retaining the conditional LAW not just the mean); (b) entropic OT for
##   d>1; (c) boundary split-guard at NIE^W=0 (WASSERSTEIN-INTEGRATION.md / BOUNDARY-SPEC.md).

if (!requireNamespace("S7", quietly = TRUE)) stop("S7 required")

## ======================= transport + EIF engine (VERIFIED, d=1) =====================
## 1-D quadratic Wasserstein via quantile transport (Brenier monotone rearrangement).
.w2_1d <- function(y1, y2, ngrid = 4096L) {
  u <- (seq_len(ngrid) - 0.5) / ngrid
  sqrt(mean((stats::quantile(y1, u, names = FALSE, type = 7) -
             stats::quantile(y2, u, names = FALSE, type = 7))^2))
}

## Kantorovich-potential influence function of W2^2(nu_a, nu_b) in its FIRST argument,
## evaluated at corner-a draws `ya` (transport toward corner-b draws `yb`).
## phi(x) = x^2 - 2 * int_0^x T(t) dt, T = Qb o Fa; IF = phi(x) - E_a[phi]. (Verified FACT5.)
.potIF <- function(ya, yb, ngrid = 4096L) {
  u  <- (seq_len(ngrid) - 0.5) / ngrid
  Qa <- stats::quantile(ya, u, names = FALSE, type = 7)
  Qb <- stats::quantile(yb, u, names = FALSE, type = 7)
  Tmap <- stats::approxfun(Qa, Qb, rule = 2)                 # monotone optimal map a->b
  grid <- seq(min(ya), max(ya), length.out = 4096L)
  Tg   <- Tmap(grid); dt <- diff(grid)
  cumI <- c(0, cumsum((utils::head(Tg, -1) + utils::tail(Tg, -1)) / 2 * dt))  # int from grid[1]
  Ifun <- stats::approxfun(grid, cumI, rule = 2)
  I0   <- Ifun(0)
  p <- ya^2 - 2 * (Ifun(ya) - I0)                            # phi(ya)
  p - mean(p)
}

## Point estimate + EIF-based SE from corner pseudo-samples (3 independent samples, d=1).
## Returns NIE^W, NDE^W, TE^W, synergy, pmedW, and (if eif=TRUE) the EIF SE of pmedW.
.pmedW_engine <- function(nu11, nu10, nu00, eif = TRUE) {
  NIE <- .w2_1d(nu11, nu10); NDE <- .w2_1d(nu10, nu00); TEW <- .w2_1d(nu11, nu00)
  D <- NIE + NDE; P <- NIE / D
  out <- list(NIE_W = NIE, NDE_W = NDE, TE_W = TEW, synergy = NIE + NDE - TEW, pmedW = P, se = NA_real_)
  if (eif) {
    n <- min(length(nu11), length(nu10), length(nu00))
    # phi_P contributions per corner sample (phi_{W2}=potIF/(2 W2); phi_P=((1-P)phi_NIE-P phi_NDE)/D)
    c11 <- ((1 - P) * .potIF(nu11, nu10) / (2 * NIE)) / D
    c10 <- ((1 - P) * .potIF(nu10, nu11) / (2 * NIE) - P * .potIF(nu10, nu00) / (2 * NDE)) / D
    c00 <- (-P * .potIF(nu00, nu10) / (2 * NDE)) / D
    out$se <- sqrt((stats::var(c11) + stats::var(c10) + stats::var(c00)) / n)
  }
  out
}

## ===================== corner-law estimator: MC g-computation (d=1) ==================
## Parametric g-computation of the four corner laws nu(a,a') for a continuous 1-D outcome.
## Reuses the family's nuisance structure (M-model, Y-model) but draws the conditional LAW.
## (Cross-fit + triply-robust corner-law EIF is the Claude-Code promotion step; see header.)
.corner_laws_gcomp <- function(d, covars, n_mc = 20000L, binY = FALSE) {
  cf  <- paste(covars, collapse = " + ")
  m_M <- stats::lm(stats::as.formula(paste("M ~ A +", cf)), data = d)
  f_Y <- stats::as.formula(paste("Y ~ A * M +", cf))
  m_Y <- stats::glm(f_Y, data = d, family = if (binY) stats::binomial() else stats::gaussian())
  sM  <- stats::sigma(m_M); sY <- if (binY) NA_real_ else stats::sigma(m_Y)
  Cdraw <- d[sample.int(nrow(d), n_mc, replace = TRUE), covars, drop = FALSE]
  corner <- function(a, ap) {
    nd_M <- data.frame(A = ap, Cdraw); muM <- stats::predict(m_M, newdata = nd_M)
    Mstar <- muM + stats::rnorm(n_mc, 0, sM)
    nd_Y <- data.frame(A = a, M = Mstar, Cdraw); muY <- stats::predict(m_Y, newdata = nd_Y, type = "response")
    if (binY) stats::rbinom(n_mc, 1, muY) else muY + stats::rnorm(n_mc, 0, sY)
  }
  list(nu11 = corner(1, 1), nu10 = corner(1, 0), nu00 = corner(0, 0))
}

## =============================== S7 result class ====================================
WassersteinPmedResult <- S7::new_class(
  "WassersteinPmedResult", package = "probmed",
  properties = list(
    pmedW = S7::class_numeric, pmedW_ci = S7::class_numeric, pmedW_se = S7::class_numeric,
    NIE_W = S7::class_numeric, NDE_W = S7::class_numeric, TE_W = S7::class_numeric,
    synergy = S7::class_numeric, method = S7::class_character,
    n = S7::class_integer, ci_level = S7::class_numeric,
    call = S7::new_property(class = S7::class_any, default = NULL)
  ),
  validator = function(self) {
    if (length(self@NIE_W) && length(self@NDE_W) && (self@NIE_W + self@NDE_W) < 1e-8)
      warning("NIE^W + NDE^W near 0: P_med^W unstable (no-mediation boundary); run the boundary guard.")
  }
)

## =============================== generic + method ===================================
wasserstein_pmed <- S7::new_generic("wasserstein_pmed", dispatch_args = "object")

## d: data.frame with columns A (0/1), M, Y, and `covars`.
## method: "eif" (transport-layer EIF SE on the g-comp corner samples) or "bootstrap".
S7::method(wasserstein_pmed, S7::class_data.frame) <-
function(object, covars = "C", n_mc = 20000L, binY = FALSE,
         method = c("eif", "bootstrap"), n_boot = 200L, ci_level = 0.95, seed = 1L, ...) {
  method <- match.arg(method); set.seed(seed)
  cl <- .corner_laws_gcomp(object, covars, n_mc, binY)
  eng <- .pmedW_engine(cl$nu11, cl$nu10, cl$nu00, eif = (method == "eif"))
  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  if (method == "bootstrap") {
    bs <- vapply(seq_len(n_boot), function(b) {
      db <- object[sample.int(nrow(object), replace = TRUE), , drop = FALSE]
      clb <- .corner_laws_gcomp(db, covars, n_mc, binY)
      .pmedW_engine(clb$nu11, clb$nu10, clb$nu00, eif = FALSE)$pmedW
    }, numeric(1))
    se <- stats::sd(bs, na.rm = TRUE)
    ci <- stats::quantile(bs, c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2), names = FALSE, na.rm = TRUE)
  } else {
    se <- eng$se; ci <- c(eng$pmedW - z * se, eng$pmedW + z * se)
  }
  WassersteinPmedResult(
    pmedW = eng$pmedW, pmedW_ci = pmin(pmax(ci, 0), 1), pmedW_se = se,
    NIE_W = eng$NIE_W, NDE_W = eng$NDE_W, TE_W = eng$TE_W, synergy = eng$synergy,
    method = method, n = nrow(object), ci_level = ci_level, call = match.call()
  )
}

if (sys.nframe() == 0L) {
  set.seed(1); n <- 3000
  C <- stats::rnorm(n); A <- stats::rbinom(n, 1, plogis(0.3 * C))
  M <- 0.4 + 0.7 * A + 0.5 * C + stats::rnorm(n)
  Y <- 0.2 + 0.6 * A + 0.8 * M + 0.3 * C + stats::rnorm(n)         # additive (no AxM): reduction case
  d <- data.frame(A, M, Y, C)
  r <- wasserstein_pmed(d, covars = "C", n_mc = 20000L, method = "eif")
  cat(sprintf("smoke g-comp: P_med^W=%.3f  se=%.3f  NIE^W=%.3f NDE^W=%.3f synergy=%.3f\n",
              r@pmedW, r@pmedW_se, r@NIE_W, r@NDE_W, r@synergy))
}
