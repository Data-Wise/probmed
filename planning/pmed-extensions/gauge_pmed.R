## gauge_pmed.R — general cross-fit one-step gauge-calibrated proportion mediated.
## Phase 3: (a) general density ratio via classifier; (b) S7 result class + generic.
## Drop-in candidate for probmed/R/ (matches probmed S7 conventions). Not yet wired
## into NAMESPACE/DESCRIPTION — promote deliberately (probmed is CRAN-blocked on medfit).
##
## Estimand: corner means theta(a,a') = E[Y(a, M(a'))]; OE=t11-t00, IDE=t10-t00,
##   IIE=t01-t00, R=OE-IDE-IIE (A x M interaction), P_med=IIE/OE, W=R/OE.
## Triply-robust EIF (Tchetgen Tchetgen & Shpitser 2012); inference by ratio identity.
##
## General nuisances (no Gaussian-M assumption):
##   pi_a(C)=P(A=a|C); q_a(M,C)=P(A=a|M,C);  mu(a,m,c)=E[Y|a,m,c]
##   density ratio f(M|a',C)/f(M|a,C) = [q_{a'}/q_a] * [pi_a/pi_{a'}]   (Bayes)
##   eta_{a,a'}(C)=E[mu(a,M,C)|A=a',C] via nested regression on C within {A=a'}.

if (!requireNamespace("S7", quietly = TRUE)) stop("S7 required")

## ---- result class (probmed S7 convention) ----
GaugePmedResult <- S7::new_class(
  "GaugePmedResult", package = "probmed",
  properties = list(
    p_med = S7::class_numeric, p_med_ci = S7::class_numeric,
    W = S7::class_numeric, W_ci = S7::class_numeric,
    W_se = S7::class_numeric, W_p = S7::class_numeric,
    OE = S7::class_numeric, IDE = S7::class_numeric,
    IIE = S7::class_numeric, R = S7::class_numeric,
    theta = S7::class_numeric, method = S7::class_character,
    n = S7::class_integer, ci_level = S7::class_numeric,
    call = S7::new_property(class = S7::class_any, default = NULL)
  ),
  validator = function(self) {
    if (length(self@OE) && abs(self@OE) < 1e-8)
      warning("OE near 0: W=R/OE is unstable; report unnormalized R.")
  }
)

## ---- engine (cross-fit one-step); covars = character vector of covariate names ----
.gauge_fit <- function(d, K, binY, covars) {
  cf <- paste(covars, collapse = " + ")
  f_pi <- stats::as.formula(paste("A ~", cf))
  f_q  <- stats::as.formula(paste("A ~ M +", cf))
  f_om <- stats::as.formula(paste("Y ~ A * M +", cf))
  n <- nrow(d); folds <- sample(rep(1:K, length.out = n))
  nm <- c("11","10","01","00"); cor <- list(c(1,1),c(1,0),c(0,1),c(0,0))
  phi <- matrix(0, n, 4, dimnames = list(NULL, nm))
  for (k in 1:K) {
    tr <- d[folds != k, , drop = FALSE]; te_i <- which(folds == k); te <- d[te_i, , drop = FALSE]
    pim <- stats::glm(f_pi, data = tr, family = stats::binomial())
    qm  <- stats::glm(f_q,  data = tr, family = stats::binomial())
    om  <- stats::glm(f_om, data = tr,
                      family = if (binY) stats::binomial() else stats::gaussian())
    p1 <- stats::predict(pim, newdata = te, type = "response"); pa <- function(z) ifelse(z==1, p1, 1-p1)
    q1 <- stats::predict(qm,  newdata = te, type = "response"); qa <- function(z) ifelse(z==1, q1, 1-q1)
    mu <- function(z) stats::predict(om, newdata = transform(te, A = z), type = "response")
    for (j in 1:4) {
      a <- cor[[j]][1]; ap <- cor[[j]][2]
      muAM <- mu(a)
      ratio <- (qa(ap)/qa(a)) * (pa(a)/pa(ap))            # f(M|ap,C)/f(M|a,C)
      ## eta via nested regression of mu(a,M,C) on covars within {A=ap}, fit on TRAIN
      muAM_tr <- stats::predict(om, newdata = transform(tr, A = a), type = "response")
      sub <- tr$A == ap
      etam <- stats::lm(stats::reformulate(covars, "yy"),
                        data = cbind(data.frame(yy = muAM_tr[sub]), tr[sub, covars, drop = FALSE]))
      eta <- stats::predict(etam, newdata = te)
      phi[te_i, j] <- (te$A==a)/pa(a) * ratio * (te$Y - muAM) +
                      (te$A==ap)/pa(ap) * (muAM - eta) + eta
    }
  }
  phi
}

## ---- generic + method ----
ward_residual <- S7::new_generic("ward_residual", dispatch_args = "object")

S7::method(ward_residual, S7::class_data.frame) <- function(object, covars = "C", K = 5L,
                                                            ci_level = 0.95, seed = 1L, ...) {
  stopifnot(all(c("A","M","Y") %in% names(object)), all(covars %in% names(object)))
  set.seed(seed)
  binY <- all(object$Y %in% 0:1)
  phi <- .gauge_fit(object, K, binY, covars); n <- nrow(object)
  th <- colMeans(phi)
  OE <- th["11"]-th["00"]; IDE <- th["10"]-th["00"]; IIE <- th["01"]-th["00"]
  R <- OE-IDE-IIE; Pmed <- IIE/OE; W <- R/OE
  pOE<-phi[,"11"]-phi[,"00"]; pIDE<-phi[,"10"]-phi[,"00"]; pIIE<-phi[,"01"]-phi[,"00"]; pR<-pOE-pIDE-pIIE
  se <- function(x) stats::sd(x)/sqrt(n); zc <- stats::qnorm(1-(1-ci_level)/2)
  seP <- se((pIIE - Pmed*pOE)/OE); seW <- se((pR - W*pOE)/OE); z <- W/seW
  GaugePmedResult(
    p_med = unname(Pmed), p_med_ci = unname(c(Pmed-zc*seP, Pmed+zc*seP)),
    W = unname(W), W_ci = unname(c(W-zc*seW, W+zc*seW)), W_se = unname(seW),
    W_p = unname(2*stats::pnorm(-abs(z))),
    OE = unname(OE), IDE = unname(IDE), IIE = unname(IIE), R = unname(R),
    theta = th, method = "onestep-crossfit", n = as.integer(n),
    ci_level = ci_level, call = match.call()
  )
}

## ---- print ----
S7::method(print, GaugePmedResult) <- function(x, ...) {
  cat("Gauge-calibrated proportion mediated (", x@method, ", n=", x@n, ")\n", sep="")
  cat(sprintf("  P_med = %.3f  [%.3f, %.3f]\n", x@p_med, x@p_med_ci[1], x@p_med_ci[2]))
  cat(sprintf("  W=R/OE = %.3f  [%.3f, %.3f]  (p=%.3g)\n", x@W, x@W_ci[1], x@W_ci[2], x@W_p))
  cat(sprintf("  OE=%.3f  IDE=%.3f  IIE=%.3f  R=%.3f\n", x@OE, x@IDE, x@IIE, x@R))
  if (abs(x@W) > 0.1) cat("  ! |W| large: additive split unreliable; interpret P_med with care.\n")
  invisible(x)
}
