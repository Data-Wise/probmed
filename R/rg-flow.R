#' Renormalization-Group Flow of the Proportion Mediated
#'
#' @description
#' S7 class for the **scale-indexed proportion mediated**
#' `g_l = NIE_l / TE_l` estimated on data coarse-grained to aggregation scale
#' `l` (block size `m`). Under the coarse-graining operator `T_m` -- which
#' replaces each within-cluster block of `m` units by its mean -- `g_l` flows
#' with `m`, and the flow is governed by the mediation beta-function.
#'
#' This is the classical **proportion mediated** (a ratio of effects), *not*
#' the package's `P_med` (a probability in `[0, 1]`); the two are different
#' quantities and are deliberately named apart. See `vignette("probmed")` for
#' `P_med`.
#'
#' @details
#' `g_l` is a smooth functional of five second moments of `(A, M, Y)`, so the
#' efficient influence function is available in closed form and the plug-in
#' estimator is already efficient in the linear-Gaussian case (no cross-fitting
#' required). The sampling unit is the **cluster**, not the unit: influence
#' contributions are aggregated to the cluster before their variance is taken.
#'
#' @param g Numeric: estimated proportion mediated at this scale.
#' @param g_ci Numeric length-2: cluster-robust EIF Wald confidence interval.
#' @param g_se Numeric: cluster-robust standard error.
#' @param a_path,b_path,TE Numeric: the `a`, `b`, and total-effect components.
#' @param grad_norm Numeric: Euclidean norm of the gradient of `g` in the five
#'   moments; bounds the bias induced by a moment perturbation of a given size.
#' @param phi_cluster Numeric: per-cluster influence contributions, in cluster
#'   order. Retained so two scales estimated on the same clusters can be
#'   contrasted by [rg_flow_contrast()].
#' @param m Numeric: aggregation scale (block size); `1` is the unit scale.
#' @param n_units Integer: coarse-grained units retained at this scale.
#' @param n_clusters Integer: clusters contributing at this scale.
#' @param ci_level Numeric: confidence level.
#' @param call Call: original call.
#'
#' @export
RgFlowResult <- S7::new_class(
  "RgFlowResult", package = "probmed",
  properties = list(
    g           = S7::class_numeric,
    g_ci        = S7::class_numeric,
    g_se        = S7::class_numeric,
    a_path      = S7::class_numeric,
    b_path      = S7::class_numeric,
    TE          = S7::class_numeric,
    grad_norm   = S7::class_numeric,
    phi_cluster = S7::class_numeric,
    m           = S7::class_numeric,
    n_units     = S7::class_integer,
    n_clusters  = S7::class_integer,
    ci_level    = S7::class_numeric,
    call        = S7::new_property(class = S7::class_any, default = NULL)
  ),
  validator = function(self) {
    if (length(self@TE) && abs(self@TE) < 1e-8) {
      return("TE is numerically zero: g = NIE/TE is not a well-posed ratio here.")
    }
    NULL
  }
)

## ---------------------------------------------------------------------------
## internals
## ---------------------------------------------------------------------------

#' Estimand and its gradient in the five second moments
#'
#' `g = (a * b) / TE` with `a = sAM/sAA`, `b = (sMY*sAA - sAY*sAM)/D`,
#' `TE = sAY/sAA`, `D = sAA*sMM - sAM^2`. `sigma_YY` drops out. The five
#' partials are computed symbolically (CAS-verified) rather than numerically.
#'
#' @param Z Numeric matrix with columns `A`, `M`, `Y`.
#' @return List with the moments, path components, `g`, and its gradient.
#' @noRd
.rg_parts <- function(Z) {
  mu <- colMeans(Z)
  S <- stats::cov(Z)
  sAA <- S[1, 1]; sAM <- S[1, 2]; sAY <- S[1, 3]; sMM <- S[2, 2]; sMY <- S[2, 3]
  D <- sAA * sMM - sAM^2
  grad <- c(
    AA =  (sAM^2 * (sAY * sMM - sAM * sMY)) / (sAY * D^2),
    AM =  (sAA * (sAA * sMM * sMY + sAM^2 * sMY - 2 * sAM * sAY * sMM)) / (sAY * D^2),
    AY = -(sAA * sAM * sMY) / (sAY^2 * D),
    MM = -(sAA * sAM * (sAA * sMY - sAM * sAY)) / (sAY * D^2),
    MY =  (sAA * sAM) / (sAY * D)
  )
  a <- sAM / sAA
  b <- (sMY * sAA - sAY * sAM) / D
  TE <- sAY / sAA
  list(
    mu = mu, sig = c(AA = sAA, AM = sAM, AY = sAY, MM = sMM, MY = sMY),
    a = a, b = b, TE = TE, g = (a * b) / TE,
    grad = grad, grad_norm = sqrt(sum(grad^2))
  )
}

#' Unit-level influence contributions
#' @noRd
.rg_phi_units <- function(Z, parts) {
  cen <- sweep(Z, 2, parts$mu)
  cp <- cbind(
    AA = cen[, 1]^2, AM = cen[, 1] * cen[, 2], AY = cen[, 1] * cen[, 3],
    MM = cen[, 2]^2, MY = cen[, 2] * cen[, 3]
  )
  as.vector(sweep(cp, 2, parts$sig) %*% parts$grad)
}

#' Aggregate unit influences to the cluster (the i.i.d. sampling unit)
#' @noRd
.rg_phi_cluster <- function(phi, cluster) {
  cluster <- droplevels(as.factor(cluster))
  s <- tapply(phi, cluster, sum)
  s / (length(phi) / length(s))
}

## ---------------------------------------------------------------------------
## exported API
## ---------------------------------------------------------------------------

#' Coarse-grain multilevel data to aggregation scale `m`
#'
#' @description
#' Applies the coarse-graining operator `T_m`: within each cluster, units are
#' randomly partitioned into blocks of size `m` and each block is replaced by
#' its mean. The cluster structure is preserved, so every scale keeps the same
#' clusters and two scales remain comparable cluster-by-cluster.
#'
#' @details
#' Blocks are complete by construction: within a cluster of size `n_s`, the
#' final incomplete block is discarded, so the number of coarse-grained units
#' is `sum(floor(n_s / m))`, not `n / m`. Clusters smaller than `m` drop out
#' entirely at that scale. The composition property
#' `T_{m2}(T_{m1}(.)) = T_{m1 * m2}(.)` holds exactly on the units a partition
#' retains.
#'
#' This function calls [sample.int()] and is therefore random. It never sets a
#' seed; seed the caller (or wrap in `withr::with_seed()`) for reproducibility.
#'
#' @param data Data frame with the treatment, mediator, outcome, and cluster
#'   columns.
#' @param m Positive integer block size. `m = 1` returns the unit scale
#'   unchanged.
#' @param treatment,mediator,outcome Character: column names of `A`, `M`, `Y`.
#' @param cluster Character: column name of the cluster identifier.
#'
#' @return List with `Z` (numeric matrix of coarse-grained units, columns
#'   `A`, `M`, `Y`) and `cluster` (their cluster labels).
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(
#'   A = rnorm(200), M = rnorm(200), Y = rnorm(200),
#'   sch = rep(1:20, each = 10)
#' )
#' cg <- rg_coarse_grain(d, m = 5, cluster = "sch")
#' nrow(cg$Z)  # 40 = 20 clusters * floor(10 / 5)
#'
#' @export
rg_coarse_grain <- function(data, m, treatment = "A", mediator = "M",
                            outcome = "Y", cluster = "cluster") {
  stopifnot(is.data.frame(data), length(m) == 1L, m >= 1L, m == as.integer(m))
  cols <- c(treatment, mediator, outcome)
  missing_cols <- setdiff(c(cols, cluster), names(data))
  if (length(missing_cols)) {
    stop("column(s) not found in `data`: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }
  m <- as.integer(m)
  if (m == 1L) {
    Z <- as.matrix(data[, cols, drop = FALSE])
    colnames(Z) <- c("A", "M", "Y")
    return(list(Z = Z, cluster = data[[cluster]]))
  }
  parts <- lapply(split(data, data[[cluster]]), function(s) {
    n <- nrow(s)
    if (n < m) return(NULL)
    idx <- sample.int(n)
    nb <- n %/% m
    s2 <- s[idx[seq_len(nb * m)], cols, drop = FALSE]
    blk <- rep(seq_len(nb), each = m)
    means <- rowsum(as.matrix(s2), blk) / m
    list(Z = means, cluster = rep(s[[cluster]][1L], nb))
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) {
    stop("no cluster has at least `m` = ", m, " units.", call. = FALSE)
  }
  Z <- do.call(rbind, lapply(parts, `[[`, "Z"))
  colnames(Z) <- c("A", "M", "Y")
  list(Z = Z, cluster = unlist(lapply(parts, `[[`, "cluster"), use.names = FALSE))
}

#' Efficient-influence-function estimator of the scale-indexed proportion mediated
#'
#' @description
#' Estimates `g_l = NIE_l / TE_l` on coarse-grained data by plug-in, with a
#' **cluster-robust** standard error from the efficient influence function.
#' The plug-in is already efficient in the linear-Gaussian case, so no
#' cross-fitting is performed.
#'
#' @param Z Numeric matrix with columns `A`, `M`, `Y` (coarse-grained units, as
#'   returned by [rg_coarse_grain()]).
#' @param cluster Cluster labels for the rows of `Z`. The i.i.d. sampling unit
#'   is the cluster, not the row; passing unit-level labels would understate
#'   every standard error.
#' @param m Numeric: the aggregation scale `Z` was produced at (recorded, not
#'   recomputed).
#' @param ci_level Numeric: confidence level for the Wald interval.
#'
#' @return An [RgFlowResult].
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(
#'   A = rnorm(400), M = rnorm(400), Y = rnorm(400),
#'   sch = rep(1:40, each = 10)
#' )
#' cg <- rg_coarse_grain(d, m = 2, cluster = "sch")
#' fit <- g_hat_eif(cg$Z, cg$cluster, m = 2)
#' fit@g
#'
#' @export
g_hat_eif <- function(Z, cluster, m = 1, ci_level = 0.95) {
  Z <- as.matrix(Z)
  if (ncol(Z) != 3L) stop("`Z` must have exactly three columns (A, M, Y).", call. = FALSE)
  if (nrow(Z) != length(cluster)) {
    stop("`cluster` must have one label per row of `Z`.", call. = FALSE)
  }
  if (length(ci_level) != 1L || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single value in (0, 1).", call. = FALSE)
  }
  parts <- .rg_parts(Z)
  pj <- .rg_phi_cluster(.rg_phi_units(Z, parts), cluster)
  se <- sqrt(mean(pj^2) / length(pj))
  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  RgFlowResult(
    g = parts$g, g_ci = c(parts$g - z * se, parts$g + z * se), g_se = se,
    a_path = parts$a, b_path = parts$b, TE = parts$TE,
    # names(pj) are the cluster labels; rg_flow_contrast() matches on them, so
    # they must survive into the result (do NOT strip with as.vector()).
    grad_norm = parts$grad_norm, phi_cluster = pj,
    m = m, n_units = nrow(Z), n_clusters = length(pj),
    ci_level = ci_level, call = match.call()
  )
}

#' Paired contrast of the proportion mediated between two scales
#'
#' @description
#' Contrasts `g` at two aggregation scales estimated from the **same** clusters
#' by differencing their influence functions cluster-by-cluster. Because the
#' shared sampling noise cancels, the paired standard error is materially
#' smaller than the naive one that treats the two scales as independent.
#'
#' @details
#' The gain is largest for two nearby scales, where the influence functions are
#' most similar cluster-by-cluster; the endpoint contrast (unit scale versus
#' fully aggregated) is the worst case for pairing.
#'
#' @param fit1,fit2 [RgFlowResult] objects estimated on the same clusters. The
#'   contrast is `fit2@g - fit1@g`. Clusters present in only one of the two are
#'   dropped.
#' @param ci_level Numeric: confidence level.
#'
#' @return List with `delta`, `se_paired`, `se_naive`, `gain`
#'   (`se_naive / se_paired`), `z`, `p_value`, `ci`, and `n_clusters`.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(
#'   A = rnorm(400), M = rnorm(400), Y = rnorm(400),
#'   sch = rep(1:40, each = 10)
#' )
#' f1 <- g_hat_eif(rg_coarse_grain(d, 1, cluster = "sch")$Z,
#'                 rg_coarse_grain(d, 1, cluster = "sch")$cluster, m = 1)
#' cg2 <- rg_coarse_grain(d, 2, cluster = "sch")
#' f2 <- g_hat_eif(cg2$Z, cg2$cluster, m = 2)
#' rg_flow_contrast(f1, f2)$se_paired
#'
#' @export
rg_flow_contrast <- function(fit1, fit2, ci_level = 0.95) {
  if (!S7::S7_inherits(fit1, RgFlowResult) || !S7::S7_inherits(fit2, RgFlowResult)) {
    stop("`fit1` and `fit2` must both be RgFlowResult objects.", call. = FALSE)
  }
  p1 <- fit1@phi_cluster
  p2 <- fit2@phi_cluster
  if (is.null(names(p1)) || is.null(names(p2))) {
    if (length(p1) != length(p2)) {
      stop("influence functions have different lengths and no cluster names to ",
           "match on; re-estimate both scales on the same clusters.", call. = FALSE)
    }
    common <- seq_along(p1)
  } else {
    common <- intersect(names(p1), names(p2))
    if (!length(common)) {
      stop("the two fits share no clusters.", call. = FALSE)
    }
  }
  dphi <- p2[common] - p1[common]
  J <- length(dphi)
  delta <- fit2@g - fit1@g
  se_paired <- sqrt(mean(dphi^2) / J)
  se_naive <- sqrt(fit1@g_se^2 + fit2@g_se^2)
  z <- delta / se_paired
  zc <- stats::qnorm(1 - (1 - ci_level) / 2)
  list(
    delta = delta, se_paired = se_paired, se_naive = se_naive,
    gain = se_naive / se_paired, z = z,
    p_value = 2 * stats::pnorm(-abs(z)),
    ci = c(delta - zc * se_paired, delta + zc * se_paired),
    n_clusters = J
  )
}
