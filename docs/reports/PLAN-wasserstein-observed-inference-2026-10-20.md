# Observed-data Wasserstein inference repair

Prespecified before implementation/validation, 2026-10-20. Authorized by the
author in response to the PR #54 inference blocker. No P_med, I6/I7, release,
or submission gate changes. Work stays in PR #55 before dev -> main.

## Scope and derivation

The legacy `eif` API name will denote a model-based observed-data influence
function, not a claim of nonparametric efficiency. Assume iid observed rows,
correct Gaussian homoscedastic mediator regression, and a correct Gaussian
homoscedastic or Bernoulli-logit outcome regression with A*M and additive
covariates. C is standardized over its empirical distribution, not fixed.
Require finite moments, nonsingular regression information, positive residual
variances, and both Wasserstein components bounded away from zero for Wald
inference. The existing identification assumptions are not weakened.

Fit the two regressions and use maximum-likelihood residual variances (mean
squared residual, not degrees-of-freedom corrected sigma squared). Stack their
row-level score influences, retaining within-row covariance. For Gaussian
outcomes each conditional corner is normal with mean
mu_Y(a,0,C) + b(a)*mu_M(a',C) and variance v_Y + b(a)^2*v_M.
Integrate the normal-mixture quantile distances deterministically. For each
quantile Q(u), implicit differentiation gives

    IF_Q,i = grad_theta Q * IF_theta,i - (F(Q|C_i) - u)/f(Q).

The second term is the empirical covariate-distribution contribution. Integrate
2*(Q_j-Q_k)*(IF_Qj-IF_Qk) to obtain IF for squared distance, then apply the square
root and ratio chain rules jointly. For Bernoulli outcomes integrate the
logistic-normal probability, use W2^2 = abs(p_j-p_k), and apply the same chain.
SE is sd(row influence)/sqrt(n_observed). No simulated draw count enters it.
Bootstrap must refit this same deterministic estimator, rather than add corner
simulation noise. Retain n_mc as an ignored compatibility argument and document
the numerical/API change. Both near-zero components trigger a diagnostic
bootstrap fallback; this does not establish boundary coverage validity.

Primary methodological support: Saul & Hudgens (2020), JSS 92(2),
doi:10.18637/jss.v092.i02, PDF pp.1-3, Section 2, equations (1)-(2), read directly.
It supports stacked estimating-equation sandwich inference, not this new
Wasserstein-specific derivation or efficiency. The derivation above is ours
and requires numerical and repeated-sampling validation.

## Test-first and independent checks

1. Regression: fixed observed rows, n_mc 2,000 versus 200,000 must give identical
   point estimates and SE; old public implementation must fail this test.
2. Differentiate weighted-refit point estimates numerically and compare with
   analytic row influences (relative tolerance 0.01, absolute 0.005).
3. Check Gaussian closed-form and Bernoulli distance identities; quadrature
   doubling changes P by <1e-4 and SE by <1% on selected regular fixtures.
4. Compare row bootstrap SE (at least 400 replicates) with IF SE within 20% on
   selected regular Gaussian and binary fixtures. This is not a coverage proof.
5. Check invalid data, singular fits, missingness, and both boundary directions.
   Mutation: omit model or covariate uncertainty and require a targeted failure.

## ADEMP repeated-sampling validation

Aim: diagnose bias, SE calibration, and Wald coverage in regular supported
models. No claim of uniform validity or extension to DR/multidimensional paths.
Cells: Gaussian and Bernoulli outcomes; n=500 and 1,000; 2,000 independent
datasets per cell. Seed 20261020, independent deterministic per-replicate seeds.
C~N(0,1); A|C~Bernoulli(expit(0.3*C));
M=0.7*A+0.5*C+epsilon_M, epsilon_M~N(0,1).
Gaussian Y=0.6*A+0.8*M+0.2*A*M+0.3*C+epsilon_Y, epsilon_Y~N(0,1).
Binary Y~Bernoulli(expit(-0.5+0.6*A+0.8*M+0.2*A*M+0.3*C)).
Gaussian truth uses the independent normal-marginal formula; binary truth
uses independent adaptive normal integration. Target P=W11,10/(W11,10+W10,00).

Save all replicate estimates, SE, coverage and failures; report elapsed time,
R/package versions, bias MCSE, empirical SD, mean estimated SE, coverage and
sqrt(p*(1-p)/2000) coverage MCSE (~0.0049 at p=.95). Prespecified engineering
diagnostics: zero failed fits; absolute bias/empirical SD <0.10; mean SE/SD in
[0.90,1.10]; coverage in [0.92,0.98]. Investigate any violation, without moving
these tolerances after observing results. These are scoped repair diagnostics,
not a replacement for the original exact-Gaussian gate or a publication gate.
Full package local CI and actual pkgdown inspection remain required before
merge; no status/protection switch until success on the exact PR head.
