# Supported-model guard and stochastic controls (repair tasks T6, T7).

guard_data <- function(n = 120, seed = 11) {
  set.seed(seed)
  X <- stats::rbinom(n, 1, 0.5)
  C <- stats::rnorm(n)
  M <- 0.5 * X + 0.3 * C + stats::rnorm(n)
  Y <- 0.2 * X + 0.4 * M + 0.3 * C + stats::rnorm(n)
  data.frame(X = X, M = M, Y = Y, C = C,
             Ycount = stats::rpois(n, exp(0.2 + 0.1 * M)),
             Ybin = stats::rbinom(n, 1, stats::plogis(-0.2 + 0.5 * M)))
}

expect_rng_untouched <- function(expr) {
  set.seed(99)
  before <- .Random.seed
  force(expr)
  expect_identical(.Random.seed, before)
}

test_that("unsupported families are rejected before any random draw", {
  d <- guard_data()
  expect_rng_untouched(expect_error(
    pmed(Ycount ~ X + M + C, formula_m = M ~ X + C, data = d,
         treatment = "X", mediator = "M", family_y = stats::poisson(),
         method = "plugin"),
    "poisson\\(log\\) outcome"
  ))
  expect_error(
    pmed(Y ~ X + M + C, formula_m = M ~ X + C, data = transform(d, Y = Y + 10),
         treatment = "X", mediator = "M",
         family_y = stats::gaussian(link = "log"), method = "plugin"),
    "gaussian\\(log\\) outcome"
  )
  expect_error(
    pmed(Ybin ~ X + M + C, formula_m = M ~ X + C, data = d,
         treatment = "X", mediator = "M",
         family_y = stats::binomial(link = "cloglog"), method = "plugin"),
    "binomial\\(cloglog\\) outcome"
  )
})

test_that("the guard also covers extracted models (MediationData dispatch)", {
  d <- guard_data()
  ex <- extract_mediation(
    stats::glm(M ~ X + C, data = d), treatment = "X", mediator = "M",
    model_y = stats::glm(Ycount ~ X + M + C, data = d, family = stats::poisson()),
    data = d
  )
  expect_rng_untouched(expect_error(pmed(ex, method = "parametric_bootstrap",
                                         n_boot = 5), "poisson"))
})

test_that("supported links still run", {
  d <- guard_data()
  for (link in c("logit", "probit")) {
    r <- pmed(Ybin ~ X + M, formula_m = M ~ X, data = d,
              treatment = "X", mediator = "M",
              family_y = stats::binomial(link = link), method = "plugin",
              n_sim = 2000, seed = 1)
    expect_true(r@estimate >= 0 && r@estimate <= 1)
  }
})

test_that("interactions and transformed treatment/mediator terms are rejected", {
  d <- guard_data()
  expect_rng_untouched(expect_error(
    pmed(Y ~ X + M + C + X:M, formula_m = M ~ X + C, data = d,
         treatment = "X", mediator = "M", method = "plugin"),
    "Unsupported term\\(s\\): `X:M`"
  ))
  expect_error(
    pmed(Y ~ X + M + I(M^2) + C, formula_m = M ~ X + C, data = d,
         treatment = "X", mediator = "M", method = "plugin"),
    "`I\\(M\\^2\\)`"
  )
  expect_error(
    pmed(Y ~ X + M + C, formula_m = M ~ log(X + 1) + C, data = d,
         treatment = "X", mediator = "M", method = "plugin"),
    "`log\\(X \\+ 1\\)`"
  )
  # Transformed covariates are not the treatment or mediator: allowed.
  r <- pmed(Y ~ X + M + log(abs(C) + 1), formula_m = M ~ X + log(abs(C) + 1),
            data = d, treatment = "X", mediator = "M", method = "plugin",
            n_sim = 2000, seed = 1)
  expect_true(is.finite(r@estimate))
})

test_that("binary models with covariates are rejected until standardization exists", {
  d <- guard_data()
  expect_rng_untouched(expect_error(
    pmed(Ybin ~ X + M + C, formula_m = M ~ X + C, data = d,
         treatment = "X", mediator = "M", family_y = stats::binomial(),
         method = "plugin"),
    "covariates \\(`C`\\) with a binary"
  ))
})

test_that("offsets, missing rows, non-numeric mediators and non-0/1 responses are rejected", {
  d <- guard_data()
  expect_rng_untouched(expect_error(
    pmed(Y ~ X + M + offset(C), formula_m = M ~ X, data = d,
         treatment = "X", mediator = "M", method = "plugin"),
    "offset"
  ))
  d_na <- d
  d_na$Y[3] <- NA
  expect_rng_untouched(expect_error(
    pmed(Y ~ X + M + C, formula_m = M ~ X + C, data = d_na,
         treatment = "X", mediator = "M", method = "nonparametric_bootstrap",
         n_boot = 5),
    "missing values.*`Y`"
  ))
  d_prop <- d
  d_prop$Ybin <- d_prop$Ybin * 0.5 + 0.25
  expect_error(suppressWarnings(
    pmed(Ybin ~ X + M, formula_m = M ~ X, data = d_prop, treatment = "X",
         mediator = "M", family_y = stats::binomial(), method = "plugin")),
    "coded 0/1"
  )
})

test_that("rank-deficient designs are rejected", {
  d <- guard_data()
  d$C2 <- 2 * d$C
  expect_error(
    pmed(Y ~ X + M + C + C2, formula_m = M ~ X + C, data = d,
         treatment = "X", mediator = "M", method = "plugin"),
    "aliased.*`y_C2`"
  )
})

test_that("invalid controls fail early with informative errors", {
  d <- guard_data()
  call_pmed <- function(...) pmed(Y ~ X + M, formula_m = M ~ X, data = d,
                                  treatment = "X", mediator = "M", ...)
  expect_error(call_pmed(method = "plugin", n_sim = 0), "`n_sim`")
  expect_error(call_pmed(method = "plugin", n_sim = 10.5), "`n_sim`")
  expect_error(call_pmed(method = "plugin", ci_level = 1), "`ci_level`")
  expect_error(call_pmed(method = "parametric_bootstrap", n_boot = -1), "`n_boot`")
  expect_error(call_pmed(method = "plugin", seed = "a"), "`seed`")
  expect_error(call_pmed(method = "plugin", x_value = NA_real_), "`x_ref`")
})

test_that("n_sim reaches the core on every single-mediator stochastic path", {
  d <- guard_data()
  seen <- integer(0)
  core <- .pmed_core_simple
  testthat::local_mocked_bindings(.pmed_core_simple = function(..., n_sim) {
    seen <<- c(seen, n_sim)
    core(..., n_sim = n_sim)
  })
  for (m in c("plugin", "parametric_bootstrap", "nonparametric_bootstrap")) {
    seen <- integer(0)
    pmed(Y ~ X + M, formula_m = M ~ X, data = d, treatment = "X",
         mediator = "M", method = m, n_boot = 3, n_sim = 321, seed = 1)
    expect_true(length(seen) >= 1, label = m)
    expect_true(all(seen == 321), label = m)
  }
  # Positive control: the spy sees the default when n_sim is not supplied.
  seen <- integer(0)
  pmed(Y ~ X + M, formula_m = M ~ X, data = d, treatment = "X",
       mediator = "M", method = "parametric_bootstrap", n_boot = 2, seed = 1)
  expect_true(all(seen == 10000))
})

test_that("an explicit seed reproduces every stochastic method", {
  d <- guard_data()
  for (m in c("plugin", "parametric_bootstrap", "nonparametric_bootstrap")) {
    r1 <- pmed(Y ~ X + M, formula_m = M ~ X, data = d, treatment = "X",
               mediator = "M", method = m, n_boot = 20, n_sim = 2000, seed = 4)
    r2 <- pmed(Y ~ X + M, formula_m = M ~ X, data = d, treatment = "X",
               mediator = "M", method = m, n_boot = 20, n_sim = 2000, seed = 4)
    expect_identical(r1@estimate, r2@estimate, label = m)
    expect_identical(r1@ci_lower, r2@ci_lower, label = m)
  }
})

test_that("larger n_sim reduces Monte Carlo spread of the plug-in", {
  d <- guard_data()
  spread <- function(n_sim) {
    stats::sd(vapply(1:30, function(s) pmed(
      Y ~ X + M, formula_m = M ~ X, data = d, treatment = "X", mediator = "M",
      method = "plugin", n_sim = n_sim, seed = s
    )@estimate, numeric(1)))
  }
  expect_lt(spread(20000), spread(500))
})

test_that("the extracted-model guard catches what the formula guard cannot reach", {
  d <- guard_data()
  ex <- extract_mediation(stats::glm(M ~ X + C, data = d), treatment = "X",
                          mediator = "M", model_y = stats::glm(Y ~ X + M + C, data = d),
                          data = d)
  expect_true(.check_supported_extract(ex))
  # NULL families count as gaussian(identity).
  expect_true(.supported_family(NULL))
  expect_identical(.family_label(NULL), "gaussian(identity)")
  # Transformed mediator term arriving through an extracted model.
  bad <- ex
  bad@outcome_predictors <- c(bad@outcome_predictors, "I(M^2)")
  expect_error(.check_supported_extract(bad), "Unsupported term\\(s\\): `I\\(M\\^2\\)`")
  # Mediator missing from the outcome model.
  nomed <- ex
  nomed@outcome_predictors <- setdiff(nomed@outcome_predictors, "M")
  expect_error(.check_supported_extract(nomed), "mediator as a predictor of the outcome")
  # Non-numeric mediator in the stored data.
  chr <- ex
  dd <- d
  dd$M <- as.character(round(dd$M, 2))
  chr@data <- dd
  expect_error(.check_supported_extract(chr), "numeric `M`")
})
