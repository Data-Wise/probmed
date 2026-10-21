# Deterministic g-computation and observed-row sandwich inference. The legacy
# public method name "eif" does not imply nonparametric efficiency.
.wasserstein_rule <- function(order, normal = FALSE) {
  j <- seq_len(order - 1L)
  off <- if (normal) sqrt(j) else j / sqrt(4 * j^2 - 1)
  mat <- matrix(0, order, order)
  mat[cbind(j, j + 1L)] <- off
  mat[cbind(j + 1L, j)] <- off
  eig <- eigen(mat, symmetric = TRUE)
  list(
    x = if (normal) eig$values else (eig$values + 1) / 2,
    w = eig$vectors[1L, ]^2
  )
}

# binY matches the legacy public argument.
.wasserstein_observed <- function(d, covars = "C", binY = FALSE, # nolint: object_name_linter.
                                  influence = TRUE, order = 128L,
                                  weights = NULL) {
  n <- nrow(d)
  vars <- c("A", "M", "Y", covars)
  if (
    n < 10L || anyNA(d[, vars, drop = FALSE]) ||
      any(!is.finite(d$M)) || any(!is.finite(d$Y)) ||
      !all(d$A %in% 0:1) || length(unique(d$A)) != 2L ||
      (binY && !all(d$Y %in% 0:1))
  ) {
    stop("Wasserstein g-computation requires complete finite data and binary A (and binary Y when binY=TRUE).")
  }
  if (is.null(weights)) weights <- rep(1 / n, n)
  if (
    length(weights) != n || any(!is.finite(weights)) ||
      any(weights <= 0) || abs(sum(weights) - 1) > 1e-8
  ) {
    stop("weights must be positive probabilities, one per observed row.")
  }
  if (influence && max(abs(weights - 1 / n)) > 1e-12) {
    stop("Influence calculation requires equal observed-row weights.")
  }
  fm <- stats::reformulate(c("A", covars), "M")
  fy <- stats::reformulate(c("A * M", covars), "Y")
  # Put weights in the formula environment rather than relying on caller lookup.
  fit_m <- stats::lm(fm, d, weights = weights, na.action = stats::na.fail)
  fit_y <- if (binY) {
    stats::glm(fy, d,
      family = stats::binomial(),
      weights = weights * n, na.action = stats::na.fail
    )
  } else {
    stats::lm(fy, d, weights = weights, na.action = stats::na.fail)
  }
  bm <- stats::coef(fit_m)
  by <- stats::coef(fit_y)
  if (any(!is.finite(c(bm, by))) || (binY && !fit_y$converged)) {
    stop("Wasserstein working models must be full rank and converged.")
  }
  em <- stats::residuals(fit_m)
  ey <- d$Y - stats::fitted(fit_y)
  vm <- sum(weights * em^2)
  vy <- if (binY) 0 else sum(weights * ey^2)
  if (!is.finite(vm + vy) || vm <= 1e-12 || (!binY && vy <= 1e-12)) {
    stop("Wasserstein working models require positive residual variances.")
  }
  pm <- length(bm)
  py <- length(by)
  p <- pm + py + 1L + !binY
  im <- seq_len(pm)
  iy <- pm + seq_len(py)
  ivm <- pm + py + 1L
  ivy <- p
  if (influence) {
    xm <- stats::model.matrix(fit_m)
    xy <- stats::model.matrix(fit_y)
    wm <- crossprod(xm) / n
    wy <- crossprod(xy, xy * if (binY) {
      stats::fitted(fit_y) * (1 - stats::fitted(fit_y))
    } else {
      1
    }) / n
    psi <- matrix(0, n, p)
    psi[, im] <- (xm * em) %*% solve(wm)
    psi[, iy] <- (xy * ey) %*% solve(wy)
    psi[, ivm] <- em^2 - vm
    if (!binY) psi[, ivy] <- ey^2 - vy
  }
  design <- function(fit, nd) {
    stats::model.matrix(
      stats::delete.response(stats::terms(fit)), nd,
      contrasts.arg = fit$contrasts, xlev = fit$xlevels
    )
  }
  rule <- .wasserstein_rule(if (binY) max(40L, order %/% 2L) else order,
    normal = binY
  )
  corner <- function(a, ap) {
    nd <- d
    nd$A <- ap
    xm <- design(fit_m, nd)
    mum <- drop(xm %*% bm)
    nd$A <- a
    nd$M <- 0
    x0 <- design(fit_y, nd)
    nd$M <- 1
    slope <- design(fit_y, nd) - x0
    b <- drop(slope %*% by)
    mu <- drop(x0 %*% by) + b * mum
    dmu <- matrix(0, n, p)
    dmu[, im] <- xm * b
    dmu[, iy] <- x0 + slope * mum
    if (binY) {
      prob <- rep(0, n)
      deriv <- matrix(0, n, p)
      for (k in seq_along(rule$x)) {
        z <- rule$x[k]
        eta <- mu + b * sqrt(vm) * z
        pk <- stats::plogis(eta)
        dp <- pk * (1 - pk)
        g <- dmu
        g[, iy] <- g[, iy, drop = FALSE] + slope * sqrt(vm) * z
        g[, ivm] <- b * z / (2 * sqrt(vm))
        prob <- prob + rule$w[k] * pk
        deriv <- deriv + rule$w[k] * g * dp
      }
      avg <- sum(weights * prob)
      phi <- if (influence) drop(psi %*% colSums(deriv * weights)) + prob - avg else NULL
      return(list(q = avg, phi = phi))
    }
    variance <- vy + b^2 * vm
    s <- sqrt(variance)
    dv <- matrix(0, n, p)
    dv[, iy] <- 2 * b * vm * slope
    dv[, ivm] <- b^2
    dv[, ivy] <- 1
    qs <- numeric(length(rule$x))
    phis <- if (influence) matrix(0, n, length(qs)) else NULL
    for (k in seq_along(qs)) {
      u <- rule$x[k]
      ends <- range(mu + s * stats::qnorm(u))
      q <- if (diff(ends) < 1e-12) {
        mean(ends)
      } else {
        stats::uniroot(
          function(y) sum(weights * stats::pnorm((y - mu) / s)) - u,
          ends,
          tol = 1e-10
        )$root
      }
      qs[k] <- q
      if (influence) {
        z <- (q - mu) / s
        density <- stats::dnorm(z) / s
        f <- sum(weights * density)
        grad <- colSums((dmu * density + dv * (density * z / (2 * s))) * weights) / f
        phis[, k] <- drop(psi %*% grad) - (stats::pnorm(z) - u) / f
      }
    }
    list(q = qs, phi = phis)
  }
  corners <- list(corner(1, 1), corner(1, 0), corner(0, 0))
  distance <- function(j, k) {
    delta <- corners[[j]]$q - corners[[k]]$q
    squared <- if (binY) abs(delta) else sum(rule$w * delta^2)
    w <- sqrt(squared)
    phi <- NULL
    if (influence && w > 1e-10) {
      dp <- corners[[j]]$phi - corners[[k]]$phi
      phi <- if (binY) {
        sign(delta) * dp / (2 * w)
      } else {
        drop(dp %*% (rule$w * delta)) / w
      }
    }
    list(w = w, phi = phi)
  }
  nie <- distance(1, 2)
  nde <- distance(2, 3)
  te <- distance(1, 3)
  total <- nie$w + nde$w
  if (!is.finite(total) || total < 1e-10) {
    stop("Wasserstein proportion is undefined when both components vanish.")
  }
  prop <- nie$w / total
  phi <- if (influence && !is.null(nie$phi) && !is.null(nde$phi)) {
    ((1 - prop) * nie$phi - prop * nde$phi) / total
  } else {
    NULL
  }
  list(
    pmedW = prop, NIE_W = nie$w, NDE_W = nde$w, TE_W = te$w,
    synergy = max(0, total - te$w),
    se = if (is.null(phi)) NA_real_ else stats::sd(phi) / sqrt(n),
    influence = phi
  )
}
