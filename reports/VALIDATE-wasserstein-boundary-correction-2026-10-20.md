# Boundary-aware Wasserstein correction research

**Decision: keep PRs #54/#55 and the local-status switch blocked.** The
candidate is research code only; no production estimator or frozen diagnostic
has changed. Design was committed as 00c5380 before implementation, and the
implementation was frozen as fa353d5 before the fresh simulations.
See [prespecified design](PLAN-wasserstein-boundary-correction-2026-10-20.md).

## What was tested

Two damped second-order curvature corrections use the original shared-model
contrast covariance, with fixed damping k=4 (primary) and k=9 (sensitivity).
The correction is applied on the log-odds scale to remain bounded without
clipping. It addresses only transformation curvature, not the entire bias
of the estimated corner probabilities. Gaussian estimation is unchanged.

Nineteen targeted checks pass: an independent finite-difference Hessian,
the small-variance expansion, bounds, endpoint/undefined cases, scaling and
contrast-swap identities, and independent binary prototype agreement with
production point estimates and SEs. No symbolic proof backend was used.
The binary implementation shares quadrature nodes with production but not its
corner/IF implementation. This is not a fully independent numerical oracle;
the earlier adaptive-integration oracle tests remain separate evidence.
Replacing the correction with the uncorrected ratio in an in-memory mutation
causes the small-variance expansion test to fail (0.5227744249 versus
0.5227744172 at tolerance 1e-10). The source file was not mutated.

The study includes 2,000 fresh datasets in each of seven cells: the four
original Gaussian/binary n=500/1,000 regular cells and three binary n=500
boundary/weak-component stress cells. Seeds use 20361020 rather than the
earlier 20261020 base. All candidates use the original first-order IF SE, as
prespecified. These internal-helper simulations do not test the public
boundary-bootstrap fallback. Any change to a candidate after seeing these
results requires another independent validation set.

## Completed fresh results

All 14,000 datasets completed with zero failed fits in 269.56 seconds under
R 4.6.1. There are 34,000 saved method-by-dataset rows, because the binary
datasets are evaluated by three paired methods. The candidate source SHA-256
is `7226c6e216f21697349bc031ef29abff1a839a25166aa9e3a087b8240adc89ff`.

The primary k=4 candidate, paired with unchanged Gaussian estimation, passes
**all four original regular-cell diagnostics**. This is a scoped success, not
a claim that true bias satisfies a confidence bound or that boundary inference
is valid. The unchanged thresholds remain exactly as specified.

| Cell | n | Method | Bias | Bias/SD | SE/SD | Coverage | Frozen diagnostics |
| --- | ---: | --- | ---: | ---: | ---: | ---: | --- |
| Gaussian regular | 500 | unchanged | -0.002984 | 0.0526 | 0.9998 | 0.9465 | Pass |
| Binary regular | 500 | primary k=4 | 0.005168 | 0.0767 | 1.0368 | 0.9675 | Pass |
| Gaussian regular | 1,000 | unchanged | -0.000756 | 0.0187 | 0.9928 | 0.9505 | Pass |
| Binary regular | 1,000 | primary k=4 | 0.000207 | 0.0048 | 0.9908 | 0.9620 | Pass |
| Binary zero direct | 500 | primary k=4 | -0.327306 | 3.1588 | 3.0503 | 0.5020 | Fail |
| Binary zero indirect | 500 | primary k=4 | 0.238426 | 2.6138 | 2.6750 | 0.5605 | Fail |
| Binary weak direct | 500 | primary k=4 | -0.108120 | 1.0232 | 2.9054 | 0.7855 | Fail |

For binary n=500 regular, raw bias/SD is 0.129162 (fails); primary k=4 reduces
RMSE from 0.069558 to 0.067534. Its bias MCSE is 0.001506 and coverage MCSE
is 0.003965. At n=1,000 these MCSEs are 0.000971 and 0.004275. Full MCSEs
for every method and cell are retained in the summary CSV.

The prespecified sensitivity candidate k=9 narrowly fails binary n=500 bias:
bias/SD=0.100275 against <0.10 (bias MCSE 0.001529). Do not round this to a
pass, select a new k, or treat this fresh set as training data for another
independently confirmed candidate. No estimates are outside [0,1] for either
damped candidate, but both fail all three stress cells. The primary correction
increases RMSE relative to raw in each stress cell. Keeping estimates bounded
has not solved either boundary bias or the observed Wald-coverage failure.

Near-zero contrasts (minimum abs(contrast)/SE <2) occur in 22.65% of regular
binary n=500 datasets and 2.65% at n=1,000, versus 95.25%-96.15% of the three
stress cells. The public API's 5% relative-distance fallback criterion is a
different diagnostic; these percentages do not validate or redefine it.

Saved [replicates](data/wasserstein-boundary-2026-10-20/replicates.csv),
[summary](data/wasserstein-boundary-2026-10-20/summary.csv), and
[provenance](data/wasserstein-boundary-2026-10-20/provenance.txt) preserve the
seeds, versions, every result and failure field. Replicate file SHA-256:
`cc46f5223dee85da3373698caae82458ec019789fc7c2acbb847a4bc55206a55`.
Reproduce from the frozen source with
`Rscript scripts/validate-wasserstein-boundary-correction.R OUTPUT_DIRECTORY`.
An independent CSV audit recomputed all 17 summary rows, verified 14,000
unique new seeds with no overlap with the earlier study, and checked every
saved failure field. The archived replicate hash matches the generated file.

**Conclusion:** retain the primary candidate as evidence of improved
regular-cell performance, but do not adopt it as a general boundary repair.
Next research should construct inference from the joint probability contrasts
before the nonsmooth square-root/ratio map (for example, confidence-set
projection), rather than assuming a bounded point correction repairs Wald
coverage. That alternative is proposed, not implemented or validated here.

## Why bounded correction does not establish boundary inference

At a population zero indirect probability contrast, a regular fitted contrast
typically fluctuates at order n^(-1/2). Taking sqrt(abs(contrast)) therefore
produces an indirect distance of order n^(-1/4), not n^(-1/2). With a fixed
positive direct component, the proportion inherits that nonregular rate.
An analogous argument applies to 1-P at a zero direct component. A regular
Wald argument cannot simply be continued to these parameter values.

For positive-semidefinite contrast covariance V and k>0, the damped ratios
satisfy 0<=rj<=1/k and abs(c)<=1/(4*k), hence abs(shift)<=9/(16*k).
The latter follows from abs(V12)<=sqrt(V11*V22) and
abs(d)/(d^2+k*v)<=1/(2*sqrt(k*v)), with zero-variance cases understood by
continuity. Thus the log-odds correction is bounded and cannot remove the
n^(-1/4) boundary behavior by itself. This is an analytic explanation of the
candidate's limitation, not a uniform theorem about the full fitted estimator.

## Source support and limits

Jiao and Han, *Bias Correction with Jackknife, Bootstrap, and Taylor Series*,
[arXiv:1709.06183v4](https://arxiv.org/pdf/1709.06183), PDF introduction
(pp.1-2), Section II.A/Theorem 2 (pp.4-5), and Section II.C, equation (76)
(pp.9-10), were read directly. Their binomial-functional analysis distinguishes
pointwise expansion arguments from broader bias/variance behavior and
illustrates why generic correction is not automatically safe. Its Taylor
discussion motivates checking smoothness. It does not derive or endorse our
damped Wasserstein correction, prove its IF, or establish boundary coverage.
The initial older HTML access attempt failed; the v4 PDF was accessible.

## Integration hold

Passing the original four-cell diagnostics is necessary, not sufficient for
production adoption. Boundary-aware point stabilization must not be described
as boundary-valid inference. Corrected-estimator inference review, public-path
validation, documentation, and full exact-commit local CI remain outstanding.
No original P_med/I6/I7 gate, release hold, or submission hold is changed.
Both PRs were verified OPEN and DRAFT after the study to make the scientific
merge hold explicit. Main's existing required check and hosted workflows
remain in place. No full local software-CI pass is claimed for this research.
