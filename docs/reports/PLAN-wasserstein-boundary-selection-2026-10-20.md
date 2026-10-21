# Frozen boundary-selection investigation

Research only. Both PRs #54 and #55 remain blocked. No production adoption,
gate change, merge, or required-status switch is authorized by this experiment.

## Candidate and justification

For signed contrasts d and their observed-data covariance V, classify a contrast
as weak when abs(d_j)/sqrt(V_jj) <= lambda. If exactly one is weak, select its
axis: report P=0 for j=1, P=1 for j=2. Otherwise retain the previous k=4 damped
curvature estimate. Record both-weak cases rather than assigning an arbitrary
value at the undefined population double-zero ratio.

Primary lambda=qnorm(.975); sensitivity lambda=qnorm(.9975). These are fixed
before fresh results. The latter targets .005 null activation, not a tuned
simulation optimum. Both candidates must be reported, including failures.

Estimate uncertainty by applying the entire rule to 32768 common antithetic
Gaussian contrast draws (seed 816). For exactly one selected axis j, center
at d - V[,j]*d[j]/V[j,j], the covariance-aware constrained projection. Otherwise
center at d. Retain V. The resulting SD is a Gaussian-experiment proxy, NOT a
uniformly valid bootstrap SE. Selection discontinuities and near-zero alternatives
may invalidate its calibration. Compare with independent 262144 draws within
10% relative error on fixed fixtures before launching the study.

Invert the minimum standardized squared perpendicular contrast over the two
lines with directions (p^2, plus/minus (1-p)^2), using qchisq(.95,1).
This is Fieller-style union-null test inversion, not a 95% joint ellipse.
Under a regular joint-contrast CLT, one line is the true null, so using its
chi-square-1 cutoff gives asymptotic coverage at least .95. This does not prove
uniform validity at the double-zero point. Record point-outside-interval events;
do not enlarge intervals to hide them.

## ADEMP and unchanged diagnostics

Use the identical seven DGP cells, population truth integration, and B=2000 per
cell from validate-wasserstein-joint-contrast.R. New dataset seeds:
20561020 + cell*10000 + replicate, disjoint from all previous studies.
Compare binary raw/Wald, k4/Fieller, primary selection/Fieller, and sensitivity
selection/Fieller; preserve Gaussian raw/Wald controls. Freeze source before run.
All cells retain: zero failures, absolute bias/empirical SD < .10, mean estimated
SE/empirical SD in [.90,1.10], coverage in [.92,.98]. Degenerate empirical SD is
not a pass. Report MCSE, widths, selection rates, all raw rows and provenance.
No threshold or gate retuning after seeing results.

Independent tests: direct perpendicular-contrast statistic against geometric
inversion, symmetry/scale/axis cases, scalar versus vectorized point rule,
and independent Gaussian proxy precision. These are research checks, not full
package CI. Only a passing scientific study can advance the existing CI workflow.

## Sources and limitations

von Luxburg and Franz (2009), Statistica Sinica 19, 1095–1117, sections 2–3,
https://www3.stat.sinica.edu.tw/sstest/oldpdf/A19n312.pdf: geometric ratio-test
inversion motivation; their exact normal-mean result does not establish exact
inference for fitted nonlinear contrasts here.

Pötscher and Schneider (2010), arXiv:0806.1652v3, sections 1–3,
https://arxiv.org/pdf/0806.1652: thresholding confidence sets require attention to
moving parameters and infimal coverage. No theorem there validates this ratio
or the proposed Gaussian proxy. Fresh weak-signal failures remain disqualifying.
