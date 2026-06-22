## ward_os(): cross-fit ONE-STEP gauge-calibrated proportion mediated (Phase 2).
## Corner means theta(a,a') = E[Y(a, M(a'))] (nested means). Triply-robust EIF
## (Tchetgen Tchetgen & Shpitser 2012) per corner:
##   phi_{a,a'} = 1{A=a}/pi_a(C) * f(M|a',C)/f(M|a,C) * {Y - mu(a,M,C)}
##              + 1{A=a'}/pi_{a'}(C) * { mu(a,M,C) - eta_{a,a'}(C) }
##              + eta_{a,a'}(C),   eta_{a,a'}(C)=E[mu(a,M,C)|A=a',C].
## Functionals: OE=t11-t00, IDE=t10-t00, IIE=t01-t00, R=OE-IDE-IIE,
##   P_med=IIE/OE, W=R/OE. Inference by ratio identity phi_W=(phi_R - W phi_OE)/OE.
## Prototype: binary A; continuous M (linear-Gaussian) -> closed-form density ratio.

ward_os <- function(d, K = 5, seed = 1, mc = 120) {
  set.seed(seed); n <- nrow(d)
  binY <- all(d$Y %in% 0:1)
  folds <- sample(rep(1:K, length.out = n))
  nm <- c("11","10","01","00"); cor <- list(c(1,1),c(1,0),c(0,1),c(0,0))
  phi <- matrix(0, n, 4, dimnames = list(NULL, nm))
  for (k in 1:K) {
    tr <- folds != k; te <- which(folds == k); nte <- length(te)
    pim <- glm(A ~ C, data = d[tr,], family = binomial())
    om  <- glm(Y ~ A * M + C, data = d[tr,],
               family = if (binY) binomial() else gaussian())
    mm  <- lm(M ~ A + C, data = d[tr,]); sig <- summary(mm)$sigma
    Cte <- d$C[te]; Ate <- d$A[te]; Mte <- d$M[te]; Yte <- d$Y[te]
    p1  <- predict(pim, newdata = data.frame(C = Cte), type = "response")
    pa  <- function(a) ifelse(a == 1, p1, 1 - p1)
    meanM <- function(a) predict(mm, newdata = data.frame(A = a, C = Cte))
    mu  <- function(a, m) predict(om, newdata = data.frame(A = a, M = m, C = Cte),
                                  type = "response")
    for (j in 1:4) {
      a <- cor[[j]][1]; ap <- cor[[j]][2]
      ratio <- dnorm(Mte, meanM(ap), sig) / dnorm(Mte, meanM(a), sig)
      muAM  <- mu(a, Mte)
      ## eta_{a,ap}(C) by vectorized Monte Carlo: M* ~ N(meanM(ap), sig)
      Mbig  <- as.vector(replicate(mc, rnorm(nte, meanM(ap), sig)))
      etaM  <- predict(om, newdata = data.frame(A = a, M = Mbig, C = rep(Cte, mc)),
                       type = "response")
      eta   <- rowMeans(matrix(etaM, nrow = nte))
      phi[te, j] <- (Ate == a)/pa(a) * ratio * (Yte - muAM) +
                    (Ate == ap)/pa(ap) * (muAM - eta) + eta
    }
  }
  th <- colMeans(phi)
  OE <- th["11"]-th["00"]; IDE <- th["10"]-th["00"]
  IIE <- th["01"]-th["00"]; R <- OE-IDE-IIE; Pmed <- IIE/OE; W <- R/OE
  pOE <- phi[,"11"]-phi[,"00"]; pIDE <- phi[,"10"]-phi[,"00"]
  pIIE <- phi[,"01"]-phi[,"00"]; pR <- pOE-pIDE-pIIE
  se <- function(x) sd(x)/sqrt(n)
  seP <- se((pIIE - Pmed*pOE)/OE); seW <- se((pR - W*pOE)/OE); z <- W/seW
  ci <- function(e,s) round(c(est=unname(e), se=unname(s),
                              lo=unname(e-1.96*s), hi=unname(e+1.96*s)), 4)
  list(theta = round(th,4),
       OE=unname(round(OE,4)), IDE=unname(round(IDE,4)),
       IIE=unname(round(IIE,4)), R=unname(round(R,4)),
       P_med = ci(Pmed, seP),
       W = c(ci(W, seW), z=unname(round(z,3)), p=unname(signif(2*pnorm(-abs(z)),3))))
}
