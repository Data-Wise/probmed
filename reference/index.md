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

## Modern estimands

Cross-fitted, influence-function-based estimators beyond the core P_med:
the gauge residual (W = R/OE, with weak-identification gates),
incremental P_med, the Sobol variance-scale and Wasserstein versions,
and the scale-indexed rg-flow proportion mediated with cluster-robust
inference.

- [`ward_residual()`](https://data-wise.github.io/probmed/reference/ward_residual.md)
  : Gauge-Calibrated Proportion Mediated

- [`incr_pmed()`](https://data-wise.github.io/probmed/reference/incr_pmed.md)
  : Incremental Mediated Elasticity Curve

- [`incr_sensitivity()`](https://data-wise.github.io/probmed/reference/incr_sensitivity.md)
  : M–Y Confounding Sensitivity for the Incremental Mediated Elasticity

- [`sobol_pmed()`](https://data-wise.github.io/probmed/reference/sobol_pmed.md)
  : Sobol / Variance-Scale Proportion Mediated

- [`sobol_from_theta()`](https://data-wise.github.io/probmed/reference/sobol_from_theta.md)
  : Map corner means to the Sobol variance decomposition

- [`wasserstein_pmed()`](https://data-wise.github.io/probmed/reference/wasserstein_pmed.md)
  : Wasserstein / Transport-Scale Proportion Mediated

- [`pmedW_dr()`](https://data-wise.github.io/probmed/reference/pmedW_dr.md)
  : Doubly-Robust Wasserstein Proportion Mediated

- [`pmedW_md()`](https://data-wise.github.io/probmed/reference/pmedW_md.md)
  : Multivariate Wasserstein Proportion Mediated

- [`sinkhorn_div()`](https://data-wise.github.io/probmed/reference/sinkhorn_div.md)
  : Debiased Sinkhorn Divergence

- [`rg_flow_contrast()`](https://data-wise.github.io/probmed/reference/rg_flow_contrast.md)
  : Paired contrast of the proportion mediated between two scales

- [`rg_coarse_grain()`](https://data-wise.github.io/probmed/reference/rg_coarse_grain.md)
  :

  Coarse-grain multilevel data to aggregation scale `m`

- [`g_hat_eif()`](https://data-wise.github.io/probmed/reference/g_hat_eif.md)
  : Efficient-influence-function estimator of the scale-indexed
  proportion mediated

## Sensitivity analysis

Unmeasured mediator-outcome confounding: tipping-point analysis shared
across the estimands.

- [`pmed_sensitivity()`](https://data-wise.github.io/probmed/reference/pmed_sensitivity.md)
  : M–Y Confounding Sensitivity for a Proportion-Mediated Estimand

## S7 Class Definitions

Formal class structures for mediation analysis objects. These classes
define the object-oriented architecture of the package.

- [`PmedResult()`](https://data-wise.github.io/probmed/reference/PmedResult.md)
  : P_med Result Class
- [`GaugePmedResult()`](https://data-wise.github.io/probmed/reference/GaugePmedResult.md)
  : Gauge-Calibrated Proportion Mediated Result
- [`IncrPmedResult()`](https://data-wise.github.io/probmed/reference/IncrPmedResult.md)
  : Incremental Mediated Elasticity Result
- [`SobolPmedResult()`](https://data-wise.github.io/probmed/reference/SobolPmedResult.md)
  : Sobol / Variance-Scale Proportion Mediated Result
- [`WassersteinPmedResult()`](https://data-wise.github.io/probmed/reference/WassersteinPmedResult.md)
  : Wasserstein / Transport-Scale Proportion Mediated Result
- [`RgFlowResult()`](https://data-wise.github.io/probmed/reference/RgFlowResult.md)
  : Renormalization-Group Flow of the Proportion Mediated
- [`PmedSensitivity()`](https://data-wise.github.io/probmed/reference/PmedSensitivity.md)
  : P_med Sensitivity Result (M–Y Confounding)

## Methods and Utilities

Generic methods for working with probmed objects, including print,
summary, and visualization functions.

- [`print-PmedResult`](https://data-wise.github.io/probmed/reference/print-PmedResult.md)
  : Print PmedResult
- [`summary-PmedResult`](https://data-wise.github.io/probmed/reference/summary-PmedResult.md)
  : Summary Method
- [`plot-PmedResult`](https://data-wise.github.io/probmed/reference/plot-PmedResult.md)
  : Plot Bootstrap Distribution

## Datasets

Reference data shipped with the package.

- [`multilevel_designs`](https://data-wise.github.io/probmed/reference/multilevel_designs.md)
  : Design metadata for candidate multilevel studies

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

- [`RgFlowResult()`](https://data-wise.github.io/probmed/reference/RgFlowResult.md)
  : Renormalization-Group Flow of the Proportion Mediated

- [`SobolPmedResult()`](https://data-wise.github.io/probmed/reference/SobolPmedResult.md)
  : Sobol / Variance-Scale Proportion Mediated Result

- [`WassersteinPmedResult()`](https://data-wise.github.io/probmed/reference/WassersteinPmedResult.md)
  : Wasserstein / Transport-Scale Proportion Mediated Result

- [`extract_mediation`](https://data-wise.github.io/probmed/reference/extract_mediation.md)
  : Re-export extract_mediation from medfit

- [`g_hat_eif()`](https://data-wise.github.io/probmed/reference/g_hat_eif.md)
  : Efficient-influence-function estimator of the scale-indexed
  proportion mediated

- [`incr_pmed()`](https://data-wise.github.io/probmed/reference/incr_pmed.md)
  : Incremental Mediated Elasticity Curve

- [`incr_sensitivity()`](https://data-wise.github.io/probmed/reference/incr_sensitivity.md)
  : M–Y Confounding Sensitivity for the Incremental Mediated Elasticity

- [`multilevel_designs`](https://data-wise.github.io/probmed/reference/multilevel_designs.md)
  : Design metadata for candidate multilevel studies

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

- [`pmedW_dr()`](https://data-wise.github.io/probmed/reference/pmedW_dr.md)
  : Doubly-Robust Wasserstein Proportion Mediated

- [`pmedW_md()`](https://data-wise.github.io/probmed/reference/pmedW_md.md)
  : Multivariate Wasserstein Proportion Mediated

- [`pmed_sensitivity()`](https://data-wise.github.io/probmed/reference/pmed_sensitivity.md)
  : M–Y Confounding Sensitivity for a Proportion-Mediated Estimand

- [`print-PmedResult`](https://data-wise.github.io/probmed/reference/print-PmedResult.md)
  : Print PmedResult

- [`rg_coarse_grain()`](https://data-wise.github.io/probmed/reference/rg_coarse_grain.md)
  :

  Coarse-grain multilevel data to aggregation scale `m`

- [`rg_flow_contrast()`](https://data-wise.github.io/probmed/reference/rg_flow_contrast.md)
  : Paired contrast of the proportion mediated between two scales

- [`sinkhorn_div()`](https://data-wise.github.io/probmed/reference/sinkhorn_div.md)
  : Debiased Sinkhorn Divergence

- [`sobol_from_theta()`](https://data-wise.github.io/probmed/reference/sobol_from_theta.md)
  : Map corner means to the Sobol variance decomposition

- [`sobol_pmed()`](https://data-wise.github.io/probmed/reference/sobol_pmed.md)
  : Sobol / Variance-Scale Proportion Mediated

- [`summary-PmedResult`](https://data-wise.github.io/probmed/reference/summary-PmedResult.md)
  : Summary Method

- [`ward_residual()`](https://data-wise.github.io/probmed/reference/ward_residual.md)
  : Gauge-Calibrated Proportion Mediated

- [`wasserstein_pmed()`](https://data-wise.github.io/probmed/reference/wasserstein_pmed.md)
  : Wasserstein / Transport-Scale Proportion Mediated

- [`.check_controls()`](https://data-wise.github.io/probmed/reference/dot-check_controls.md)
  : Validate the numerical controls shared by every method

- [`.check_formula_data_complete()`](https://data-wise.github.io/probmed/reference/dot-check_formula_data_complete.md)
  : Reject missing values in the user's data before any model is fitted

- [`.check_original_separation()`](https://data-wise.github.io/probmed/reference/dot-check_original_separation.md)
  : Stop if the original fit of a binary equation is separated

- [`.check_supported_data()`](https://data-wise.github.io/probmed/reference/dot-check_supported_data.md)
  : Data-level checks: missing rows, variable types, binary responses

- [`.check_supported_extract()`](https://data-wise.github.io/probmed/reference/dot-check_supported_extract.md)
  : Validate an extracted single-mediator model before computation

- [`.check_supported_formulas()`](https://data-wise.github.io/probmed/reference/dot-check_supported_formulas.md)
  : Reject structures a formula would fit but the core cannot represent

- [`.covariate_terms()`](https://data-wise.github.io/probmed/reference/dot-covariate_terms.md)
  : Covariates the model uses beyond the treatment and mediator

- [`.crn_uniforms()`](https://data-wise.github.io/probmed/reference/dot-crn_uniforms.md)
  :

  Uniforms for common random numbers in
  [`.pmed_core_simple()`](https://data-wise.github.io/probmed/reference/dot-pmed_core_simple.md)

- [`.design_column()`](https://data-wise.github.io/probmed/reference/dot-design_column.md)
  : Rebuild one model-matrix column from the stored data

- [`.design_matrices()`](https://data-wise.github.io/probmed/reference/dot-design_matrices.md)
  : Rebuild the mediator and outcome design matrices of an extracted
  model

- [`.detect_separation()`](https://data-wise.github.io/probmed/reference/dot-detect_separation.md)
  : Deterministic separation check for a binary-response design

- [`.extract_separation()`](https://data-wise.github.io/probmed/reference/dot-extract_separation.md)
  : Separation status of an extracted model's binary equations

- [`.family_label()`](https://data-wise.github.io/probmed/reference/dot-family_label.md)
  : Describe a family for an error message

- [`.is_gaussian()`](https://data-wise.github.io/probmed/reference/dot-is_gaussian.md)
  : Is a family Gaussian (NULL is treated as Gaussian)?

- [`.label_vars()`](https://data-wise.github.io/probmed/reference/dot-label_vars.md)
  : Data variables a model-matrix column label refers to

- [`.mbco_branch_lr()`](https://data-wise.github.io/probmed/reference/dot-mbco_branch_lr.md)
  : Branch likelihood-ratio statistics for H0: a = 0 and H0: b = 0

- [`.mbco_fit_pmed()`](https://data-wise.github.io/probmed/reference/dot-mbco_fit_pmed.md)
  : Signed Gaussian P_med profile with positive variance coordinates

- [`.mbco_gll()`](https://data-wise.github.io/probmed/reference/dot-mbco_gll.md)
  : Gaussian profile log-likelihood given SSE, variance V, and n rows.

- [`.mbco_invert()`](https://data-wise.github.io/probmed/reference/dot-mbco_invert.md)
  : Invert an LR test: find where excess(t) = LR(t) - crit crosses zero.

- [`.mbco_invert_verified()`](https://data-wise.github.io/probmed/reference/dot-mbco_invert_verified.md)
  : Verified bounded-resolution LR inversion for single-mediator models

- [`.mbco_ll_constrained_ie()`](https://data-wise.github.io/probmed/reference/dot-mbco_ll_constrained_ie.md)
  :

  Maximized log-likelihood under the constraint `a*b` = ie0

- [`.mbco_ll_constrained_ie_parallel()`](https://data-wise.github.io/probmed/reference/dot-mbco_ll_constrained_ie_parallel.md)
  : Maximized log-likelihood under the constraint sum(a\*b) = ie0
  (parallel).

- [`.mbco_ll_constrained_pmed()`](https://data-wise.github.io/probmed/reference/dot-mbco_ll_constrained_pmed.md)
  : Maximized log-likelihood under the constraint P_med = pnorm(qstar).

- [`.mbco_ll_constrained_pmed_parallel()`](https://data-wise.github.io/probmed/reference/dot-mbco_ll_constrained_pmed_parallel.md)
  : Maximized log-likelihood under the constraint joint P_med =
  pnorm(qstar).

- [`.mbco_ll_null()`](https://data-wise.github.io/probmed/reference/dot-mbco_ll_null.md)
  : Gaussian exact-null union likelihood

- [`.mbco_lr_at()`](https://data-wise.github.io/probmed/reference/dot-mbco_lr_at.md)
  : Profile LR statistic for H0: P_med = p (single-mediator Gaussian
  MBCO)

- [`.mbco_nlopt()`](https://data-wise.github.io/probmed/reference/dot-mbco_nlopt.md)
  :

  One nloptr run in the [`optim()`](https://rdrr.io/r/stats/optim.html)
  result shape

- [`.mbco_optim()`](https://data-wise.github.io/probmed/reference/dot-mbco_optim.md)
  : Convergence-checked nloptr minimizer with one deterministic fallback

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

- [`.mbco_quadratic()`](https://data-wise.github.io/probmed/reference/dot-mbco_quadratic.md)
  : Cached quadratic residual sum of squares for fixed nuisance designs

- [`.mbco_sse()`](https://data-wise.github.io/probmed/reference/dot-mbco_sse.md)
  : Evaluate a precomputed residual quadratic

- [`.named_or()`](https://data-wise.github.io/probmed/reference/dot-named_or.md)
  : Look up a named coefficient, returning a default when absent

- [`.needs_standardization()`](https://data-wise.github.io/probmed/reference/dot-needs_standardization.md)
  : Does this model need covariate standardization?

- [`.offending_terms()`](https://data-wise.github.io/probmed/reference/dot-offending_terms.md)
  : Model terms that involve a variable other than as its own main
  effect

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

- [`.pmed_model_boot()`](https://data-wise.github.io/probmed/reference/dot-pmed_model_boot.md)
  : Model-based response bootstrap for P_med (opt-in)

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

- [`.split_coefs()`](https://data-wise.github.io/probmed/reference/dot-split_coefs.md)
  : Split a joint coefficient vector into mediator and outcome parts

- [`.std_for()`](https://data-wise.github.io/probmed/reference/dot-std_for.md)
  : Standardization inputs for one coefficient vector, or NULLs if not
  needed

- [`.std_linear_predictors()`](https://data-wise.github.io/probmed/reference/dot-std_linear_predictors.md)
  : Per-row covariate linear predictors at x = 0 and m = 0

- [`.supported_family()`](https://data-wise.github.io/probmed/reference/dot-supported_family.md)
  : Is a (family, link) pair one the core can simulate?

- [`lm_class`](https://data-wise.github.io/probmed/reference/lm_class.md)
  : S3 Class Wrappers for S7

- [`` `%||%` ``](https://data-wise.github.io/probmed/reference/null-coalesce.md)
  : Null coalescing operator

- [`probmed-package`](https://data-wise.github.io/probmed/reference/probmed-package.md)
  [`probmed`](https://data-wise.github.io/probmed/reference/probmed-package.md)
  : probmed: Probabilistic Effect Sizes for Causal Mediation Analysis

- [`standardize`](https://data-wise.github.io/probmed/reference/standardize.md)
  : Observed-sample covariate standardization

- [`support-guard`](https://data-wise.github.io/probmed/reference/support-guard.md)
  : Supported-model guard for single-mediator P_med
