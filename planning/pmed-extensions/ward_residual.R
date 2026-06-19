## ward_residual(): gauge-calibrated proportion mediated for Proposal 4
## Returns (P_med, Ward residual W = R/OE) with bootstrap CIs, where
##   theta(a,a') = E_C E[ mu(a, M, C) ],  M ~ p(M | A=a', C)   (stochastic/interventional draw)
##   OE=theta(1,1)-theta(0,0); IDE=theta(1,0)-theta(0,0); IIE=theta(0,1)-theta(0,0)
##   R = OE - IDE - IIE  (mixed 2nd difference = direct-x-mediator interaction)
## Plug-in g-computation (Phase 1 deliverable; medrobust EIF comes in Phase 2).

## ---- core point estimator ----
.ward_point <- function(d, K = 200) {
  ## outcome + mediator working models
  om <- glm(Y ~ A * M + C, data = d,
            family = if (all(d$Y %in% 0:1)) binomial() else gaussian())
  mm <- lm(M ~ A + C, data = d); sigM <- summary(mm)$sigma
  n  <- nrow(d); C <- d$C
  ## theta(a, a'): set A=a in outcome; draw M from mediator model at A=a'
  theta <- function(a, ap) {
    mu_acc <- numeric(n)
    for (k in seq_len(K)) {
      Mstar <- predict(mm, newdata = data.frame(A = ap, C = C)) + rnorm(n, 0, sigM)
      mu_acc <- mu_acc + predict(om, newdata = data.frame(A = a, M = Mstar, C = C),
                                 type = "response")
    }
    mean(mu_acc / K)
  }
  t11 <- theta(1,1); t10 <- theta(1,0); t01 <- theta(0,1); t00 <- theta(0,0)
  OE <- t11 - t00; IDE <- t10 - t00; IIE <- t01 - t00; R <- OE - IDE - IIE
  c(OE = OE, IDE = IDE, IIE = IIE, R = R, P_med = IIE / OE, W = R / OE)
}

## ---- public function with nonparametric bootstrap ----
ward_residual <- function(d, K = 200, B = 200, seed = 1) {
  set.seed(seed)
  est <- .ward_point(d, K)
  boot <- replicate(B, .ward_point(d[sample(nrow(d), replace = TRUE), ], K))
  ci <- apply(boot, 1, quantile, c(.025, .975), na.rm = TRUE)
  list(estimate = round(est, 4),
       ci = round(t(ci), 4),
       test_W = c(W = unname(est["W"]),
                  z = unname(est["R"] / sd(boot["R", ], na.rm = TRUE)),
                  p = unname(2 * pnorm(-abs(est["R"] / sd(boot["R", ], na.rm = TRUE))))))
}
