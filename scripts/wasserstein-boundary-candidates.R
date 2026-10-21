# Research prototypes only; not sourced by the package's public API.
wb_ratio <- function(d) {
  stopifnot(length(d) == 2L, all(is.finite(d)))
  w <- sqrt(abs(d))
  if (sum(w) == 0) stop("Both zero components: ratio undefined")
  w[1] / sum(w)
}

wb_correct <- function(d, v, k = 4) {
  stopifnot(
    all(dim(v) == c(2L, 2L)), all(is.finite(v)), k > 0,
    all(diag(v) >= 0)
  )
  p <- wb_ratio(d)
  if (p == 0 || p == 1) {
    return(p)
  }
  den <- d^2 + k * diag(v)
  r <- diag(v) / den
  cross <- v[1, 2] * d[1] * d[2] / prod(den)
  shift <- (-(1 + 2 * p) * r[1] + (3 - 2 * p) * r[2]) / 8 -
    (1 - 2 * p) * cross / 4
  stats::plogis(stats::qlogis(p) - shift)
}

# Independent binary implementation, vectorizing the integration rather than
# calling the production corner/IF code. The Gaussian rule itself is shared.
wb_binary <- function(d) {
  n <- nrow(d)
  fm <- stats::lm(M ~ A + C, d)
  fy <- stats::glm(Y ~ A * M + C, d, family = stats::binomial())
  if (!fy$converged || anyNA(c(stats::coef(fm), stats::coef(fy)))) {
    stop("Working-model fit failed")
  }
  xm <- stats::model.matrix(fm)
  xy <- stats::model.matrix(fy)
  em <- stats::residuals(fm)
  py <- stats::fitted(fy)
  vm <- mean(em^2)
  if (vm <= 1e-12) stop("Mediator variance vanished")
  psi <- cbind(
    (xm * em) %*% solve(crossprod(xm) / n),
    (xy * (d$Y - py)) %*% solve(crossprod(xy, xy * py * (1 - py)) / n),
    em^2 - vm
  )
  rule <- probmed:::.wasserstein_rule(64L, normal = TRUE)
  mm <- function(fit, nd) {
    stats::model.matrix(stats::delete.response(stats::terms(fit)), nd,
      contrasts.arg = fit$contrasts, xlev = fit$xlevels
    )
  }
  corner <- function(a, ap) {
    nd <- d
    nd$A <- ap
    cm <- mm(fm, nd)
    mum <- drop(cm %*% stats::coef(fm))
    nd$A <- a
    nd$M <- 0
    x0 <- mm(fy, nd)
    nd$M <- 1
    slope <- mm(fy, nd) - x0
    b <- drop(slope %*% stats::coef(fy))
    mu <- drop(x0 %*% stats::coef(fy)) + b * mum
    eta <- outer(mu, rep(1, length(rule$x))) + outer(b * sqrt(vm), rule$x)
    pk <- stats::plogis(eta)
    prob <- drop(pk %*% rule$w)
    dp <- drop((pk * (1 - pk)) %*% rule$w)
    zdp <- drop((pk * (1 - pk)) %*% (rule$w * rule$x))
    grad <- c(
      colMeans(cm * (b * dp)),
      colMeans(x0 * dp + slope * (mum * dp + sqrt(vm) * zdp)),
      mean(b * zdp / (2 * sqrt(vm)))
    )
    list(p = mean(prob), phi = drop(psi %*% grad) + prob - mean(prob))
  }
  cs <- list(corner(1, 1), corner(1, 0), corner(0, 0))
  contrasts <- c(cs[[1]]$p - cs[[2]]$p, cs[[2]]$p - cs[[3]]$p)
  phi <- cbind(cs[[1]]$phi - cs[[2]]$phi, cs[[2]]$phi - cs[[3]]$phi)
  v <- stats::cov(phi) / n
  p <- wb_ratio(contrasts)
  grad <- p * (1 - p) * c(0.5 / contrasts[1], -0.5 / contrasts[2])
  list(
    raw = p, d = contrasts, v = v,
    se = sqrt(drop(t(grad) %*% v %*% grad)),
    min_z = min(abs(contrasts) / sqrt(diag(v)))
  )
}
