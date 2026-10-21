# Observed-sample covariate standardization

For a binary mediator or outcome, covariates do not cancel in the
cross-world contrast, so P_med is standardized over the empirical
covariate distribution of the analysis rows: each of the `n_sim` pairs
draws one row uniformly and uses its covariate values in both arms
(Definition 1 with the target held fixed at the observed rows). Gaussian
additive models are not standardized: their covariate shifts cancel
exactly.

## Details

medfit stores predictors as model-matrix column names (`"C"`, `"Fb"`,
`"log(abs(C) + 1)"`), so the design is rebuilt column by column from the
stored data.
