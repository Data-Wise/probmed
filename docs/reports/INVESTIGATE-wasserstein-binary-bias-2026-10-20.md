# Follow-up diagnosis: binary n=500 bias

The original four-cell study, source d476a0f, completed 8,000 datasets with no
failed fits. Three cells meet every prespecified engineering diagnostic. Binary
n=500 does not: bias/empirical SD is 0.142817 rather than <0.10. Its mean
SE/empirical SD is 0.988387 and coverage is 0.9630 (MCSE 0.004221). This is a
failed bias diagnostic, not a pass, despite acceptable SE and coverage results.

Before running the follow-up, specify this decomposition: replay the same 2,000
binary n=500 datasets and independently reconstruct the three fitted corner
probabilities. Check reproduction of all saved estimates within 1e-8. Let g be
the probability-to-Wasserstein-ratio map. Decompose total bias exactly as

    mean(g(p_hat)) - g(p_true)
      = [g(mean(p_hat)) - g(p_true)]
        + [mean(g(p_hat)) - g(mean(p_hat))].

Report both terms and the counts of estimated indirect/direct probability
contrasts reversing their true positive signs. This distinguishes nuisance
probability bias from curvature/near-boundary behavior of the nonlinear ratio.
It neither changes the estimator nor relaxes the original diagnostic. No main
merge or required-status switch on the strength of acceptable coverage alone.

## Completed replay and author decision

All 2,000 reconstructed estimates match the saved estimates within 1e-8.
Total bias is 0.009851444: fitted-probability bias contributes 0.001396391,
and nonlinear transformation contributes 0.008455053 (85.8%). There are zero
indirect-contrast and four direct-contrast sign reversals. This decomposition
is descriptive, not a proof of a particular bias order or a usable correction.
The author selected **keep merge blocked; investigate bias correction** on
2026-10-20. PRs #54/#55 and the required-status switch remain blocked.

## Exploratory grouped-jackknife pilot (specified before running)

Do not change the production estimator or the frozen validation thresholds.
For the first 200 previously generated binary n=500 datasets, partition rows
deterministically into 20 balanced groups. Refit after omitting each group and
compute P_corrected = 20*P_full - 19*mean(P_leave_group_out). This targets
leading smooth-estimator bias without changing the Wasserstein estimand.
Grouping is a computational pilot, not a substitute for ordinary delete-one
jackknife validation. The same deterministic row ordering is independent of the
generated iid observations; each omitted group has 25 rows.

Report paired bias, RMSE, correction variability, estimates outside [0,1],
failed fits, and nominal Wald coverage using the original IF SE solely as an
exploratory diagnostic. A first-order SE need not describe the corrected
estimator adequately in finite samples. Do not clip estimates or suppress
failures. Near-zero components make the ratio nonsmooth and can destabilize
jackknife extrapolation. This is a development sample, not independent
confirmation; 200 datasets cannot establish a small standardized-bias bound.
Any promising correction requires fresh seeds, at least 2,000 datasets per
cell, independent checks, and review of boundary handling before adoption.

## Pilot results: do not adopt unguarded jackknife correction

The 200 development datasets completed without fit errors. Raw bias was
0.003689 (MCSE 0.004799), versus corrected bias -0.002087 (MCSE 0.007936).
This small development subset does not reproduce the full 2,000-dataset bias
precisely and cannot establish a reduction. RMSE increased from 0.067798 to
0.111977; SD increased from 0.067868 to 0.112238. Mean correction was
-0.005777 (paired MCSE 0.005611); correction SD was 0.079346.
One corrected estimate was outside [0,1]. Exploratory coverage using the
uncorrected IF SE was 0.960 raw and 0.955 corrected, not evidence validating
that SE for the corrected estimator.

Replicate 166 changed from 0.868992 to 1.913700 under grouped jackknife.
A targeted ordinary delete-one replay on this same dataset gave -1.764091,
with deleted estimates ranging from 0.810100 to 0.992244. Thus replacing
grouped deletion with ordinary deletion does not remove the observed extreme
behavior. These are post-pilot diagnostics, not independent validation.

The raw IF SE for this selected dataset is 0.823595, indicating severe local
instability. The study calls the internal regular-model helper, not the public
boundary-bootstrap fallback; neither this pilot nor the original four-cell
study establishes boundary-fallback coverage. Unguarded extrapolation can
magnify nonsmooth near-zero-component behavior. Reject this candidate for
production; do not repair it by clipping or relaxing the frozen diagnostics.
Next research should explicitly distinguish regular separated components from
near-zero components before proposing a bounded, boundary-aware correction.
Full local CI remains required for any eventual production candidate, but
software CI alone cannot lift this author-selected scientific merge hold.

Scripts: `scripts/pilot-wasserstein-bias-correction.R` and
`scripts/diagnose-wasserstein-jackknife-extreme.R`. Saved pilot replicates and
summary are under `data/wasserstein-observed-2026-10-20/`.
