args <- commandArgs(TRUE)
source("inst/sim/hopper/i6/config.R")
source("inst/sim/hopper/i6/mbco_pmed.R")

mixture <- .7 * pnorm(.5 * .8 / sqrt(2 * (1 + .8^2))) + .3 * .5
slope <- .7 * .8
residual_var <- 1 + .7 * .3 * .8^2 * (1 + .5 * .5^2)
gaussian <- pnorm(.5 * slope / sqrt(2 * (slope^2 + residual_var)))
score <- function(data) {
  n <- nrow(data)
  m <- .lm.fit(cbind(1, data$X), data$M)
  y <- .lm.fit(cbind(1, data$X, data$M), data$Y)
  a <- m$coefficients[2L]
  b <- y$coefficients[3L]
  pnorm(a * b / sqrt(2 * (b^2 * sum(m$residuals^2) / (n - 2L) +
                           sum(y$residuals^2) / (n - 3L))))
}
RNGkind("L'Ecuyer-CMRG")
set.seed(20261123)
streams <- vector("list", 2000L)
streams[[1L]] <- .Random.seed
for (i in 2:2000) streams[[i]] <- parallel::nextRNGStream(streams[[i - 1L]])
one <- function(i) {
  assign(".Random.seed", streams[[i]], envir = .GlobalEnv)
  d <- SCENARIOS$s5$gen(1000L)
  mbco <- pmed_mbco_ci(d)
  estimate <- score(d)
  stopifnot(abs(mbco$estimate - estimate) < 1e-12)
  draws <- replicate(1000L, score(d[sample.int(1000L, 1000L, TRUE), ]))
  ci <- quantile(draws, c(.025, .975), names = FALSE)
  c(i = i, est = estimate, boot_sd = sd(draws), mlo = mbco$ci_lower,
    mhi = mbco$ci_upper, blo = ci[1L], bhi = ci[2L])
}
nreps <- as.integer(args[[2L]])
rows <- parallel::mclapply(seq_len(nreps), one, mc.cores = 8L)
stopifnot(all(vapply(rows, is.numeric, logical(1))))
d <- as.data.frame(do.call(rbind, rows))
stopifnot(identical(as.integer(d$i), seq_len(nreps)))
write.csv(d, args[[1L]], row.names = FALSE)
cat("targets:", mixture, gaussian, "gap:", gaussian - mixture, "\n")
print(c(bias = mean(d$est) - mixture, bias_mcse = sd(d$est) / sqrt(nreps),
        sd = sd(d$est), shift_over_sd = (gaussian - mixture) / sd(d$est),
        boot_sd = mean(d$boot_sd)))
covers <- function(target, lo, hi) is.finite(lo) & is.finite(hi) & lo <= target & target <= hi
for (prefix in c("m", "b")) {
  lo <- d[[paste0(prefix, "lo")]]; hi <- d[[paste0(prefix, "hi")]]
  x <- covers(mixture, lo, hi); y <- covers(gaussian, lo, hi)
  print(prefix)
  print(c(n = nrow(d), missing = sum(!is.finite(lo) | !is.finite(hi)),
          mixture = sum(x), gaussian = sum(y), true_only = sum(x & !y),
          working_only = sum(y & !x), mean_width = mean(hi - lo)))
  print(rbind(mixture = binom.test(sum(x), nreps)$conf.int,
              working = binom.test(sum(y), nreps)$conf.int))
}
