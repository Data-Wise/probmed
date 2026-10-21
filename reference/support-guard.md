# Supported-model guard for single-mediator P_med

The core draws Gaussian responses (identity link) or Bernoulli responses
(logit or probit link) and nothing else. Fitting a GLM is not the same
as being able to generate its response law, so every other family, link
and non-additive structure is rejected here, before any random draw.
