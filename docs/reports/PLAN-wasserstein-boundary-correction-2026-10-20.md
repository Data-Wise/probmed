# Boundary-aware correction research: frozen fresh-seed study

Author request: investigate boundary-aware bias corrections using fresh
simulations; retain the merge hold until the frozen diagnostics pass. This
research does not change production code, the estimand, or any P_med/I6 gate.

## Candidate derivation (ours, not a published recommendation)

For Bernoulli corner probabilities let d=(p11-p10,p10-p00),
P=sqrt(abs(d1))/(sqrt(abs(d1))+sqrt(abs(d2))), and V be the estimated covariance
matrix of the two estimated contrasts, including shared fitted-model and
empirical-covariate uncertainty. Away from zero contrasts,
h=logit(P)=0.5*(log(abs(d1))-log(abs(d2))). Its gradient is
a=(0.5/d1,-0.5/d2), and its Hessian is diag(-0.5/d1^2,0.5/d2^2).
The second-order transformation contribution to bias(P), divided by P*(1-P),
is 0.5*trace(((1-2*P)*a*a' + Hessian(h))*V). This omits nuisance-estimation
bias and higher-order remainder terms. It is not a complete bias expansion.

The unregularized expression can explode near zero contrasts. Define
D1=d1^2+k*V11, D2=d2^2+k*V22, r1=V11/D1, r2=V22/D2,
c=V12*d1*d2/(D1*D2), and

    shift = (-(1+2*P)*r1 + (3-2*P)*r2)/8 - (1-2*P)*c/4
    corrected = plogis(qlogis(P) - shift).

This retains the leading separated-component correction as V tends to zero,
damps singular denominators, and keeps estimates bounded without clipping.
Exact single-zero contrast estimates stay at their endpoint; both-zero
contrasts are undefined and must error. This is numerical boundary awareness,
NOT a theorem guaranteeing boundary inference or unbiasedness.

Primary candidate k=4; sensitivity candidate k=9. Freeze both now. No selection
or tuning on the confirmation results. Compare the original estimator on the
same fresh datasets. Gaussian estimation is unchanged. For all candidates use
the original first-order IF SE: the correction is O(1/n) for separated fixed
contrasts, but this SE justification is not uniform near zero. The simulation
must assess SE calibration instead of assuming it. No public bootstrap
fallback is used in these internal-helper research runs.

## Fresh validation and stress tests

Four original regular cells: Gaussian and binary, n=500/1000, 2,000 independent
datasets each. Use the previously frozen DGP, independent truth integrals and
unchanged diagnostic limits: zero failed fits, abs(bias)/empirical SD <0.10,
mean SE/SD in [0.90,1.10], coverage in [0.92,0.98]. New per-replicate seeds:
20361020 + 10000*cell + replicate. These are disjoint from the earlier study.
Report all methods, raw outcomes, bias/coverage MCSE, RMSE, correction sizes,
and near-zero-contrast frequency (abs(d)/SE(d)<2). Independent prototype
binary raw estimates and SEs must match the production helper on fixtures.

Additional binary n=500 stress cells, 2,000 datasets each: (1) zero direct
component, setting the outcome A and A:M coefficients to zero; (2) zero
indirect component, setting mediator A coefficient to zero; (3) weak direct
component, setting outcome A=0.05, A:M=0. Both-zero population components are
undefined and are tested as a mathematical input error, not a coverage cell.
Stress seeds use cell IDs 5,6,7 with the same new seed formula. These cells
must be reported even if the four original cells pass. Do not interpret
regular-cell success as validation of exact-boundary Wald inference. No merge
if a required frozen diagnostic fails; all passing is necessary, not by itself
sufficient for adopting a new production estimator or lifting other holds.

## Verification and decision

Before simulation: independent numerical Hessian check, boundedness,
contrast-swap complementarity, scaling invariance, endpoint/undefined tests,
small-V expansion, and binary prototype vs production point/SE agreement.
Record source commit and hashes. Any revision after these fresh outcomes is
a new candidate and needs another fresh simulation seed set. Do not recycle
this validation set as independent confirmation of a revised estimator.

No full software CI is claimed by this study. Production adoption would also
require corrected-estimator inference review, public-path integration,
documented limitations, and full exact-commit local CI before any merge.
