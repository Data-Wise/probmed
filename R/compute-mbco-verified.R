#' Cached quadratic residual sum of squares for fixed nuisance designs
#' @keywords internal
.mbco_quadratic <- function(design, response, predictor) {
  q <- qr(design)
  y <- qr.resid(q, response)
  x <- qr.resid(q, predictor)
  c(sum(y^2), sum(y * x), sum(x^2))
}

#' Evaluate a precomputed residual quadratic
#' @keywords internal
.mbco_sse <- function(quad, coefficient) {
  max(0, quad[1] - 2 * coefficient * quad[2] + coefficient^2 * quad[3])
}

#' Convergence-checked optimizer with one deterministic retry
#' @keywords internal
.mbco_optim <- function(start, fn) {
  first <- tryCatch(stats::optim(start, fn, method = "BFGS",
    control = list(maxit = 1500, reltol = 1e-12)), error = function(e) NULL)
  if (!is.null(first) && first$convergence == 0 && is.finite(first$value) && first$value < 1e99) return(first)
  tryCatch(stats::optim(start, fn, method = if (length(start) > 1) "Nelder-Mead" else "BFGS",
    control = list(maxit = 3000, reltol = 1e-12)), error = function(e) NULL)
}

#' Signed Gaussian P_med profile with positive variance coordinates
#' @keywords internal
.mbco_fit_pmed <- function(prep, qstar) {
  if (qstar == 0) {
    ll <- .mbco_ll_null(prep)
    return(list(ll = ll, converged = is.finite(ll), basins = list()))
  }
  target_sign <- sign(qstar * prep$delta)
  yy <- prep$quad_y[1]
  ym <- prep$quad_y[2]
  basins <- lapply(c(-1, 1), function(sa) {
    decode <- function(p) {
      Vm <- exp(p[1])
      r <- exp(p[2])
      # Profile b at fixed r = Vy/b^2. Its stationary equation is
      # n*r*b^2 + ym*b - yy = 0; each sign has one admissible root.
      disc <- sqrt(ym^2 + 4 * prep$n * r * yy)
      sb <- target_sign * sa
      b <- if (sb > 0) 2 * yy / (disc + ym) else -2 * yy / (disc - ym)
      a <- sa * abs(qstar) * sqrt(2 * (Vm + r)) / abs(prep$delta)
      list(a = a, b = b, Vm = Vm, Vy = b^2 * r)
    }
    nll <- function(p) {
      v <- decode(p)
      if (!all(is.finite(unlist(v))) || v$Vm <= 0 || v$Vy <= 0) return(1e100)
      ll <- .mbco_gll(.mbco_sse(prep$quad_m, v$a), v$Vm, prep$n) +
        .mbco_gll(.mbco_sse(prep$quad_y, v$b), v$Vy, prep$n)
      if (is.finite(ll)) -ll else 1e100
    }
    bm <- max(abs(prep$b_hat), 0.1)
    base <- log(c(prep$Vm_hat, prep$Vy_hat / bm^2))
    starts <- list(base, base + c(0, -log(4)))
    # Near the null each sign basin has two branches: a ~ 0 (b near b_hat,
    # reached from `base`) and b ~ 0 (a near a_hat). The second needs a large
    # r; holding a at |a_hat| in the constraint gives r = (a delta / q)^2 / 2 - Vm.
    r_b0 <- (prep$a_hat * prep$delta / qstar)^2 / 2 - prep$Vm_hat
    if (is.finite(r_b0) && r_b0 > 0) {
      starts[[length(starts) + 1L]] <- c(log(prep$Vm_hat), log(r_b0))
    }
    fits <- lapply(starts, function(s) .mbco_optim(s, nll))
    valid <- vapply(fits, function(o) !is.null(o) && o$convergence == 0 &&
      is.finite(o$value) && o$value < 1e99, logical(1))
    if (!any(valid)) return(list(converged = FALSE))
    fits <- fits[valid]
    best <- fits[[which.min(vapply(fits, function(o) o$value, numeric(1)))]]
    v <- decode(best$par)
    actual <- v$a * prep$delta * v$b / sqrt(2 * (v$b^2 * v$Vm + v$Vy))
    ll <- -nll(best$par)
    ok <- is.finite(actual) && abs(stats::pnorm(actual) - stats::pnorm(qstar)) <= 1e-8 &&
      ll <= prep$ll_free + 1e-6
    list(converged = ok, ll = ll, parameters = v)
  })
  ok <- all(vapply(basins, function(b) isTRUE(b$converged), logical(1)))
  list(ll = if (ok) max(vapply(basins, function(b) b$ll, numeric(1))) else NA_real_,
       converged = ok, basins = basins)
}

#' Verified bounded-resolution LR inversion for single-mediator models
#'
#' A full-domain scan at spacing no larger than step/2 detects observed
#' re-entry. Narrower components can remain undetected; this is not a proof
#' of global connectedness. Failed profiles invalidate the reported interval.
#' @keywords internal
.mbco_invert_verified <- function(est, excess, domain, step, anchors = numeric(0)) {
  failed <- function(reason = "invalid_domain_or_point", points = numeric(0)) list(bounds = c(lower = NA_real_, upper = NA_real_),
    endpoint_status = c(lower = "failed", upper = "failed"), status = "failed",
    verified = FALSE, components = matrix(numeric(0), ncol = 2),
    reason = reason, failed_profile_points = points)
  if (!all(is.finite(c(est, domain, step))) || length(domain) != 2 ||
      domain[1] >= domain[2] || step <= 0 || est < domain[1] || est > domain[2]) return(failed())
  grid <- sort(unique(c(seq(domain[1], domain[2], length.out =
    ceiling(diff(domain) / (step / 2)) + 1L), est,
    anchors[anchors >= domain[1] & anchors <= domain[2]])))
  values <- vapply(grid, function(x) tryCatch(excess(x), error = function(e) NA_real_), numeric(1))
  if (any(!is.finite(values))) return(failed("profile_failure", grid[!is.finite(values)]))
  if (values[match(est, grid)] > 0) return(failed("estimate_rejected"))
  accepted <- values <= 0
  starts <- which(accepted & !c(FALSE, utils::head(accepted, -1)))
  ends <- which(accepted & !c(utils::tail(accepted, -1), FALSE))
  components <- matrix(NA_real_, length(starts), 2, dimnames = list(NULL, c("lower", "upper")))
  statuses <- matrix("root", length(starts), 2)
  tol <- min(1e-9, step / 1000000)
  root <- function(i, j) {
    ans <- tryCatch(stats::uniroot(excess, c(grid[i], grid[j]), tol = tol), error = function(e) NULL)
    if (is.null(ans) || !is.finite(ans$f.root) || abs(ans$f.root) > 1e-5) return(NA_real_)
    # Verify signs in the original bracket and on either side of the root.
    epsilon <- min((grid[j] - grid[i]) / 4, max(10 * tol, 1e-7))
    left <- tryCatch(excess(max(grid[i], ans$root - epsilon)), error = function(e) NA_real_)
    right <- tryCatch(excess(min(grid[j], ans$root + epsilon)), error = function(e) NA_real_)
    if (!all(is.finite(c(left, right))) ||
        !((left <= 0 && right >= 0) || (left >= 0 && right <= 0))) return(NA_real_)
    ans$root
  }
  for (k in seq_along(starts)) {
    i <- starts[k]; j <- ends[k]
    if (i == 1) {
      components[k, 1] <- domain[1]; statuses[k, 1] <- "truncated"
    } else components[k, 1] <- root(i - 1, i)
    if (j == length(grid)) {
      components[k, 2] <- domain[2]; statuses[k, 2] <- "truncated"
    } else components[k, 2] <- root(j, j + 1)
  }
  if (anyNA(components)) return(failed("root_verification"))
  status <- if (nrow(components) > 1) "disconnected" else
    if (any(statuses == "truncated")) "truncated" else "ok"
  list(bounds = if (status == "ok" || status == "truncated")
    stats::setNames(as.numeric(components[1, ]), c("lower", "upper")) else c(lower = NA_real_, upper = NA_real_),
    endpoint_status = if (nrow(components) == 1) stats::setNames(statuses[1, ], c("lower", "upper")) else
      c(lower = "failed", upper = "failed"),
    status = status, verified = status == "ok", components = components,
    scan_spacing = max(diff(grid)), root_tolerance = tol)
}
