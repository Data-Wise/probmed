# Package index

## Core Functions

Primary functions for computing P_med and extracting mediation
structures from fitted statistical models.

- [`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md) :
  Compute P_med: Probabilistic Effect Size for Mediation
- [`pmed-formula`](https://data-wise.github.io/probmed/reference/pmed-formula.md)
  : Compute P_med from Formula
- [`pmed-MediationData`](https://data-wise.github.io/probmed/reference/pmed-MediationData.md)
  : Compute P_med from MediationData
- [`extract_mediation`](https://data-wise.github.io/probmed/reference/extract_mediation.md)
  : Re-export extract_mediation from medfit

## S7 Class Definitions

Formal class structures for mediation analysis objects. These classes
define the object-oriented architecture of the package.

- [`PmedResult()`](https://data-wise.github.io/probmed/reference/PmedResult.md)
  : P_med Result Class

## Methods and Utilities

Generic methods for working with probmed objects, including print,
summary, and visualization functions.

- [`print-PmedResult`](https://data-wise.github.io/probmed/reference/print-PmedResult.md)
  : Print PmedResult
- [`summary-PmedResult`](https://data-wise.github.io/probmed/reference/summary-PmedResult.md)
  : Summary Method
- [`plot-PmedResult`](https://data-wise.github.io/probmed/reference/plot-PmedResult.md)
  : Plot Bootstrap Distribution

## Internal Functions

Low-level computational functions. These functions are exported for
transparency but are not intended for direct use by end users.

- [`GaugePmedResult()`](https://data-wise.github.io/probmed/reference/GaugePmedResult.md)
  : Gauge-Calibrated Proportion Mediated Result

- [`IncrPmedResult()`](https://data-wise.github.io/probmed/reference/IncrPmedResult.md)
  : Incremental Mediated Elasticity Result

- [`PmedResult()`](https://data-wise.github.io/probmed/reference/PmedResult.md)
  : P_med Result Class

- [`PmedSensitivity()`](https://data-wise.github.io/probmed/reference/PmedSensitivity.md)
  : P_med Sensitivity Result (M–Y Confounding)

- [`SobolPmedResult()`](https://data-wise.github.io/probmed/reference/SobolPmedResult.md)
  : Sobol / Variance-Scale Proportion Mediated Result

- [`extract_mediation`](https://data-wise.github.io/probmed/reference/extract_mediation.md)
  : Re-export extract_mediation from medfit

- [`incr_pmed()`](https://data-wise.github.io/probmed/reference/incr_pmed.md)
  : Incremental Mediated Elasticity Curve

- [`plot-PmedResult`](https://data-wise.github.io/probmed/reference/plot-PmedResult.md)
  : Plot Bootstrap Distribution

- [`pmed-MediationData`](https://data-wise.github.io/probmed/reference/pmed-MediationData.md)
  : Compute P_med from MediationData

- [`pmed-ParallelMediationData`](https://data-wise.github.io/probmed/reference/pmed-ParallelMediationData.md)
  : Compute parallel (joint) P_med from a ParallelMediationData object

- [`pmed-formula`](https://data-wise.github.io/probmed/reference/pmed-formula.md)
  : Compute P_med from Formula

- [`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md) :
  Compute P_med: Probabilistic Effect Size for Mediation

- [`pmed_sensitivity()`](https://data-wise.github.io/probmed/reference/pmed_sensitivity.md)
  : M–Y Confounding Sensitivity for a Proportion-Mediated Estimand

- [`print-PmedResult`](https://data-wise.github.io/probmed/reference/print-PmedResult.md)
  : Print PmedResult

- [`sobol_from_theta()`](https://data-wise.github.io/probmed/reference/sobol_from_theta.md)
  : Map corner means to the Sobol variance decomposition

- [`sobol_pmed()`](https://data-wise.github.io/probmed/reference/sobol_pmed.md)
  : Sobol / Variance-Scale Proportion Mediated

- [`summary-PmedResult`](https://data-wise.github.io/probmed/reference/summary-PmedResult.md)
  : Summary Method

- [`ward_residual()`](https://data-wise.github.io/probmed/reference/ward_residual.md)
  : Gauge-Calibrated Proportion Mediated

- [`.is_gaussian()`](https://data-wise.github.io/probmed/reference/dot-is_gaussian.md)
  : Is a family Gaussian (NULL is treated as Gaussian)?

- [`.mbco_gll()`](https://data-wise.github.io/probmed/reference/dot-mbco_gll.md)
  : Gaussian profile log-likelihood given SSE, variance V, and n rows.

- [`.mbco_invert()`](https://data-wise.github.io/probmed/reference/dot-mbco_invert.md)
  : Invert an LR test: find where excess(t) = LR(t) - crit crosses zero.

- [`.mbco_ll_constrained_ie()`](https://data-wise.github.io/probmed/reference/dot-mbco_ll_constrained_ie.md)
  :

  Maximized log-likelihood under the constraint `a*b` = ie0.

- [`.mbco_ll_constrained_ie_parallel()`](https://data-wise.github.io/probmed/reference/dot-mbco_ll_constrained_ie_parallel.md)
  : Maximized log-likelihood under the constraint sum(a\*b) = ie0
  (parallel).

- [`.mbco_ll_constrained_pmed()`](https://data-wise.github.io/probmed/reference/dot-mbco_ll_constrained_pmed.md)
  : Maximized log-likelihood under the constraint P_med = pnorm(qstar).

- [`.mbco_ll_constrained_pmed_parallel()`](https://data-wise.github.io/probmed/reference/dot-mbco_ll_constrained_pmed_parallel.md)
  : Maximized log-likelihood under the constraint joint P_med =
  pnorm(qstar).

- [`.mbco_parallel_ll()`](https://data-wise.github.io/probmed/reference/dot-mbco_parallel_ll.md)
  : Joint Gaussian log-likelihood for the parallel free model at given
  params.

- [`.mbco_parallel_optim_best()`](https://data-wise.github.io/probmed/reference/dot-mbco_parallel_optim_best.md)
  :

  Run Nelder-Mead from several starts, returning the best maximized
  log-likelihood (`-min value`), or NA if every start fails / is
  penalized.

- [`.mbco_parallel_starts()`](https://data-wise.github.io/probmed/reference/dot-mbco_parallel_starts.md)
  : Deterministic perturbed warm starts for the parallel MBCO
  optimizers.

- [`.mbco_prep()`](https://data-wise.github.io/probmed/reference/dot-mbco_prep.md)
  : Build design matrices and free-fit summaries for the MBCO
  optimization.

- [`.mbco_prep_parallel()`](https://data-wise.github.io/probmed/reference/dot-mbco_prep_parallel.md)
  : Build design matrices and free-fit summaries for the parallel MBCO
  fit.

- [`.named_or()`](https://data-wise.github.io/probmed/reference/dot-named_or.md)
  : Look up a named coefficient, returning a default when absent

- [`.ols_sse()`](https://data-wise.github.io/probmed/reference/dot-ols_sse.md)
  : Residual sum of squares of an OLS fit of y on design matrix D.

- [`.pmed_compute()`](https://data-wise.github.io/probmed/reference/dot-pmed_compute.md)
  : Core P_med Computation Dispatcher

- [`.pmed_compute_parallel()`](https://data-wise.github.io/probmed/reference/dot-pmed_compute_parallel.md)
  : Dispatch parallel (joint) P_med by inference method.

- [`.pmed_core_simple()`](https://data-wise.github.io/probmed/reference/dot-pmed_core_simple.md)
  : Simple P_med Core Computation

- [`.pmed_mbco()`](https://data-wise.github.io/probmed/reference/dot-pmed_mbco.md)
  : MBCO Confidence Interval for P_med (Gaussian single-mediator model)

- [`.pmed_nonparametric_boot()`](https://data-wise.github.io/probmed/reference/dot-pmed_nonparametric_boot.md)
  : Nonparametric Bootstrap for P_med

- [`.pmed_parallel_boot_result()`](https://data-wise.github.io/probmed/reference/dot-pmed_parallel_boot_result.md)
  :

  Assemble a bootstrap `PmedResult` (shared by both parallel
  bootstraps).

- [`.pmed_parallel_closed()`](https://data-wise.github.io/probmed/reference/dot-pmed_parallel_closed.md)
  : Closed-form joint P_med for parallel mediators (Gaussian).

- [`.pmed_parallel_coef_names()`](https://data-wise.github.io/probmed/reference/dot-pmed_parallel_coef_names.md)
  :

  Names of the structural coefficients (a_j, b_j) in the parallel
  estimates vector / vcov. Uses the structural aliases `a{j}` / `b{j}`,
  which BOTH the lm/glm and lavaan medfit extractors expose (with full
  vcov rows) – unlike the source-specific `m{j}_<tx>` / `y_<mediator_j>`
  names, which exist only in the lm/glm extractor and would break the
  parametric bootstrap on a lavaan extract.

- [`.pmed_parallel_mbco()`](https://data-wise.github.io/probmed/reference/dot-pmed_parallel_mbco.md)
  : MBCO confidence interval for parallel joint P_med (Gaussian).

- [`.pmed_parallel_nonparametric_boot()`](https://data-wise.github.io/probmed/reference/dot-pmed_parallel_nonparametric_boot.md)
  : Nonparametric bootstrap for parallel joint P_med.

- [`.pmed_parallel_parametric_boot()`](https://data-wise.github.io/probmed/reference/dot-pmed_parallel_parametric_boot.md)
  : Parametric bootstrap for parallel joint P_med.

- [`.pmed_parallel_plugin()`](https://data-wise.github.io/probmed/reference/dot-pmed_parallel_plugin.md)
  : Plugin (point estimate only) for parallel joint P_med.

- [`.pmed_parallel_require_gaussian()`](https://data-wise.github.io/probmed/reference/dot-pmed_parallel_require_gaussian.md)
  : Require a Gaussian parallel extract (finite residual SDs present).

- [`.pmed_parametric_boot()`](https://data-wise.github.io/probmed/reference/dot-pmed_parametric_boot.md)
  : Parametric Bootstrap for P_med

- [`.pmed_plugin()`](https://data-wise.github.io/probmed/reference/dot-pmed_plugin.md)
  : Plugin Estimator (Point Estimate Only)

- [`.qr_sse()`](https://data-wise.github.io/probmed/reference/dot-qr_sse.md)
  : Residual sum of squares against a precomputed QR decomposition.

- [`lm_class`](https://data-wise.github.io/probmed/reference/lm_class.md)
  : S3 Class Wrappers for S7

- [`` `%||%` ``](https://data-wise.github.io/probmed/reference/null-coalesce.md)
  : Null coalescing operator

- [`probmed-package`](https://data-wise.github.io/probmed/reference/probmed-package.md)
  [`probmed`](https://data-wise.github.io/probmed/reference/probmed-package.md)
  : probmed: Probabilistic Effect Sizes for Causal Mediation Analysis
