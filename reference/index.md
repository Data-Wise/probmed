# Package index

## Core Functions

Primary functions for computing P_med and extracting mediation
structures from fitted statistical models.

- [`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md) :
  Compute P_med
- [`pmed-formula`](https://data-wise.github.io/probmed/reference/pmed-formula.md)
  : Compute P_med from Formula
- [`pmed-MediationExtract`](https://data-wise.github.io/probmed/reference/pmed-MediationExtract.md)
  : Compute P_med from MediationExtract
- [`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md)
  : Extract Mediation Structure
- [`extract_mediation-lm`](https://data-wise.github.io/probmed/reference/extract_mediation-lm.md)
  : Extract from GLM/LM

## S7 Class Definitions

Formal class structures for mediation analysis objects. These classes
define the object-oriented architecture of the package.

- [`MediationExtract()`](https://data-wise.github.io/probmed/reference/MediationExtract.md)
  : Mediation Extract Base Class
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

- [`MediationExtract()`](https://data-wise.github.io/probmed/reference/MediationExtract.md)
  : Mediation Extract Base Class
- [`PmedResult()`](https://data-wise.github.io/probmed/reference/PmedResult.md)
  : P_med Result Class
- [`extract_mediation-lm`](https://data-wise.github.io/probmed/reference/extract_mediation-lm.md)
  : Extract from GLM/LM
- [`extract_mediation()`](https://data-wise.github.io/probmed/reference/extract_mediation.md)
  : Extract Mediation Structure
- [`plot-PmedResult`](https://data-wise.github.io/probmed/reference/plot-PmedResult.md)
  : Plot Bootstrap Distribution
- [`pmed-MediationExtract`](https://data-wise.github.io/probmed/reference/pmed-MediationExtract.md)
  : Compute P_med from MediationExtract
- [`pmed-formula`](https://data-wise.github.io/probmed/reference/pmed-formula.md)
  : Compute P_med from Formula
- [`pmed()`](https://data-wise.github.io/probmed/reference/pmed.md) :
  Compute P_med
- [`print-PmedResult`](https://data-wise.github.io/probmed/reference/print-PmedResult.md)
  : Print PmedResult
- [`summary-PmedResult`](https://data-wise.github.io/probmed/reference/summary-PmedResult.md)
  : Summary Method
- [`.pmed_compute()`](https://data-wise.github.io/probmed/reference/dot-pmed_compute.md)
  : Core P_med Computation Dispatcher
- [`.pmed_core_simple()`](https://data-wise.github.io/probmed/reference/dot-pmed_core_simple.md)
  : Simple P_med Core Computation
- [`.pmed_nonparametric_boot()`](https://data-wise.github.io/probmed/reference/dot-pmed_nonparametric_boot.md)
  : Nonparametric Bootstrap for P_med
- [`.pmed_parametric_boot()`](https://data-wise.github.io/probmed/reference/dot-pmed_parametric_boot.md)
  : Parametric Bootstrap for P_med
- [`.pmed_plugin()`](https://data-wise.github.io/probmed/reference/dot-pmed_plugin.md)
  : Plugin Estimator (Point Estimate Only)
- [`lm_class`](https://data-wise.github.io/probmed/reference/lm_class.md)
  : S3 Class Wrappers for S7
- [`` `%||%` ``](https://data-wise.github.io/probmed/reference/null-coalesce.md)
  : Null coalescing operator
- [`probmed-package`](https://data-wise.github.io/probmed/reference/probmed-package.md)
  [`probmed`](https://data-wise.github.io/probmed/reference/probmed-package.md)
  : probmed: Probabilistic Effect Sizes for Causal Mediation Analysis
