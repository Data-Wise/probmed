# Covariance-extraction validation (repair task T8)

Date: 2026-10-14. Test: `tests/testthat/test-covariance-extraction.R`.
medfit tested: **0.7.0** (the installed version; DESCRIPTION requires >= 0.3.0).
The 0.4.0 and 0.6.0 routes named in the 2026-10-09 plan are no longer installed
here; the contract test below fails loudly if a future medfit changes the
structure, which is the property that matters for probmed.

## What medfit returns

`extract_mediation()@vcov` is block diagonal: the mediator block equals
`vcov(glm_m)`, the outcome block equals `vcov(glm_y)`, and every cross entry is
exactly 0. The aliases `a`, `b`, `c_prime` duplicate `m_<X>`, `y_<M>`, `y_<X>`
exactly, so the matrix is singular but positive semidefinite (smallest
eigenvalue ~ -5e-18); `MASS::mvrnorm` handles this, and the parametric bootstrap
selects coefficients by the `m_` / `y_` names.

Gaussian blocks use the residual-df variance (`summary(lm)$sigma^2`), not ML
SSE/n. MBCO uses ML variances (batch 2); the parametric bootstrap keeps the
glm convention. Both are deliberate.

## Disposition: zero cross block justified, model-based

The sequential likelihood f(M | X, C) f(Y | X, M, C) factorizes, so the Fisher
information, and hence the model-based covariance, has no cross block. This
holds for Gaussian and Bernoulli outcomes and for shared or differing covariate
sets. It is **not** a universal independence claim:

| Check (independent stacked sandwich, n = 20,000) | corr(â, b̂) |
|---|---|
| Gaussian Y, shared covariate | −0.014 |
| Bernoulli Y, shared covariate | 0.016 |
| Gaussian Y, outcome-only covariate | −0.010 |
| Misspecified outcome (omits (e_m² − 1)·X) — positive control | 0.424 |

Under correct specification the sandwich cross-correlation is O(n^-1/2) noise.
Under misspecification it is not, and the model-based zero block understates
the joint uncertainty. The opt-in response bootstrap (tasks T12–T13) refits
both equations on the same rows and so keeps coefficient pairing; it still
assumes the specified sequential model.

A first attempt used a finite paired-refit bootstrap correlation with a
4/√B band. It flagged 0.154 in one configuration at n = 300 because the
bootstrap correlation is itself a dataset-specific quantity with O(n^-1/2)
spread, so the band was wrong, not the extraction. The spec warns against
exactly this test.
