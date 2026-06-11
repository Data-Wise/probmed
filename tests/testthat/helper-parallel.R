# Helpers for parallel (joint) P_med tests:
#   - generate_parallel_data(): k-mediator Gaussian DGM
#   - extract_parallel(): fit + medfit::extract_mediation(structure = "parallel")
#   - oracle_parallel_pmed(): independent brute-force joint-estimand simulator
#   - oracle_parallel_mbco_ci(): full-dimensional (un-profiled) MBCO interval

#' Generate parallel-mediator Gaussian data: X -> Mj (independent), Y ~ X + sum(Mj).
#'
#' @return data.frame with columns X, M1..Mk, Y; attr "true_effects".
generate_parallel_data <- function(n = 1500,
                                   a_vec = c(0.5, 0.4),
                                   b_vec = c(0.5, 0.6),
                                   c_prime = 0.3,
                                   sigma_m = NULL,
                                   sigma_y = 1,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  k <- length(a_vec)
  stopifnot(length(b_vec) == k)
  if (is.null(sigma_m)) sigma_m <- rep(1, k)
  if (length(sigma_m) == 1L) sigma_m <- rep(sigma_m, k)

  X <- stats::rnorm(n)
  M <- matrix(NA_real_, nrow = n, ncol = k)
  for (j in seq_len(k)) {
    M[, j] <- a_vec[j] * X + stats::rnorm(n, sd = sigma_m[j])
  }
  Y <- c_prime * X + as.numeric(M %*% b_vec) + stats::rnorm(n, sd = sigma_y)

  data <- data.frame(X = X)
  for (j in seq_len(k)) data[[paste0("M", j)]] <- M[, j]
  data$Y <- Y
  attr(data, "true_effects") <- list(
    a_vec = a_vec, b_vec = b_vec, c_prime = c_prime,
    sigma_m = sigma_m, sigma_y = sigma_y,
    indirect_effect = sum(a_vec * b_vec)
  )
  data
}

#' Fit the k mediator models + outcome model and extract a ParallelMediationData.
extract_parallel <- function(data, meds = NULL) {
  if (is.null(meds)) meds <- grep("^M[0-9]+$", names(data), value = TRUE)
  k <- length(meds)
  m_models <- lapply(meds, function(m) {
    stats::lm(stats::reformulate("X", response = m), data = data)
  })
  form_y <- stats::reformulate(c("X", meds), response = "Y")
  model_y <- stats::lm(form_y, data = data)
  medfit::extract_mediation(
    object = m_models[[1]],
    model_y = model_y,
    treatment = "X",
    mediator = meds,
    mediator_models = m_models[-1],
    structure = "parallel",
    data = data
  )
}

#' Independent brute-force simulator of the joint parallel P_med estimand.
oracle_parallel_pmed <- function(a_vec, b_vec, sigma_m_vec, sigma_y, delta,
                                 n_sim = 4e6, seed = 1) {
  set.seed(seed)
  k <- length(a_vec)
  diff <- numeric(n_sim)
  for (j in seq_len(k)) {
    m_t <- stats::rnorm(n_sim, mean = a_vec[j] * delta, sd = sigma_m_vec[j])
    m_c <- stats::rnorm(n_sim, mean = 0, sd = sigma_m_vec[j])
    diff <- diff + b_vec[j] * (m_t - m_c)
  }
  yt <- diff + stats::rnorm(n_sim, 0, sigma_y)
  yc <- stats::rnorm(n_sim, 0, sigma_y)
  # equivalent to comparing two outcomes sharing the mediator-difference shift
  mean((yt) > (yc)) + 0.5 * mean(yt == yc)
}

# --- Full-dimensional MBCO oracle (no profiling), k mediators, no covariates ---

.oracle_par_ll_free <- function(data, meds) {
  ll <- 0
  for (m in meds) {
    ll <- ll + as.numeric(stats::logLik(
      stats::lm(stats::reformulate("X", response = m), data)
    ))
  }
  ll + as.numeric(stats::logLik(
    stats::lm(stats::reformulate(c("X", meds), response = "Y"), data)
  ))
}

.oracle_par_ll_constrained <- function(data, qstar, delta, meds) {
  n <- nrow(data)
  k <- length(meds)
  X <- data$X
  Mmat <- as.matrix(data[, meds, drop = FALSE])
  Y <- data$Y

  if (abs(qstar) < 1e-8) {
    ll <- 0
    for (m in meds) {
      ll <- ll + as.numeric(stats::logLik(
        stats::lm(stats::reformulate("X", response = m), data)
      ))
    }
    return(ll + as.numeric(stats::logLik(stats::lm(Y ~ X, data))))
  }

  s_target <- sign(qstar) * sign(delta)
  # Warm start from free fits.
  fm <- lapply(meds, function(m) stats::lm(stats::reformulate("X", response = m), data))
  fy <- stats::lm(stats::reformulate(c("X", meds), response = "Y"), data)
  a0 <- vapply(fm, function(f) unname(stats::coef(f)["X"]), numeric(1))
  im0 <- vapply(fm, function(f) unname(stats::coef(f)[1]), numeric(1))
  lvm0 <- vapply(fm, function(f) log(stats::sigma(f)^2), numeric(1))
  cy <- stats::coef(fy)
  iy0 <- unname(cy[1])
  cp0 <- unname(cy["X"])
  b0 <- vapply(meds, function(m) unname(cy[m]), numeric(1))
  start <- c(im0, a0, lvm0, iy0, cp0, b0)

  nll <- function(p) {
    i <- 0
    im <- p[i + seq_len(k)]; i <- i + k
    a <- p[i + seq_len(k)]; i <- i + k
    Vm <- exp(p[i + seq_len(k)]); i <- i + k
    iy <- p[i + 1L]; i <- i + 1L
    cp <- p[i + 1L]; i <- i + 1L
    b <- p[i + seq_len(k)]
    S <- sum(a * b)
    if (abs(S) > 1e-10 && sign(S) != s_target) return(1e10)
    Vy <- (delta * S)^2 / (2 * qstar^2) - sum(b^2 * Vm)
    if (!is.finite(Vy) || Vy <= 1e-8) return(1e10)
    ll <- 0
    for (j in seq_len(k)) {
      ll <- ll + sum(stats::dnorm(Mmat[, j], im[j] + a[j] * X, sqrt(Vm[j]), log = TRUE))
    }
    mu_y <- iy + cp * X + as.numeric(Mmat %*% b)
    ll <- ll + sum(stats::dnorm(Y, mu_y, sqrt(Vy), log = TRUE))
    -ll
  }
  o <- try(stats::optim(start, nll, method = "Nelder-Mead",
    control = list(maxit = 8000, reltol = 1e-11)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= 1e9) return(NA_real_)
  -o$value
}

oracle_parallel_mbco_ci <- function(data, meds, delta = 1, level = 0.95) {
  fm <- lapply(meds, function(m) stats::lm(stats::reformulate("X", response = m), data))
  fy <- stats::lm(stats::reformulate(c("X", meds), response = "Y"), data)
  a <- vapply(fm, function(f) unname(stats::coef(f)["X"]), numeric(1))
  b <- vapply(meds, function(m) unname(stats::coef(fy)[m]), numeric(1))
  Vm <- vapply(fm, function(f) stats::sigma(f)^2, numeric(1))
  Vy <- stats::sigma(fy)^2
  S <- sum(a * b)
  est <- as.numeric(stats::pnorm(delta * S / sqrt(2 * sum(b^2 * Vm) + 2 * Vy)))
  ll1 <- .oracle_par_ll_free(data, meds)
  crit <- stats::qchisq(level, df = 1)
  excess <- function(ps) {
    ps <- min(max(ps, 1e-5), 1 - 1e-5)
    ll0 <- .oracle_par_ll_constrained(data, stats::qnorm(ps), delta, meds)
    if (is.na(ll0)) return(crit + 1e3)
    -2 * (ll0 - ll1) - crit
  }
  endpoint <- function(dir) {
    inner <- est
    for (kk in seq_len(200)) {
      outer <- min(max(est + dir * 0.01 * kk, 1e-4), 1 - 1e-4)
      if (excess(outer) > 0) {
        return(tryCatch(stats::uniroot(excess, sort(c(inner, outer)), tol = 1e-4)$root,
          error = function(e) NA_real_))
      }
      inner <- outer
      if (outer <= 1e-4 || outer >= 1 - 1e-4) return(outer)
    }
    outer
  }
  c(estimate = est, lower = endpoint(-1), upper = endpoint(+1))
}
