# Frozen joint-contrast projection investigation

Research only. Both PRs #54/#55 remain drafts. No production inference change,
merge, protection switch, or hosted-workflow change is authorized by this study.

## Procedure frozen before implementation and fresh results

For binary outcomes let d=(p11-p10,p10-p00), with observed-row covariance V
from the independent `wb_binary()` prototype. Use the joint Wald region
E={(x-d)' V^-1 (x-d) <= qchisq(.95,2)}. Return the closed image of E excluding
the origin under g(x)=sqrt(abs(x1))/(sqrt(abs(x1))+sqrt(abs(x2))).
No chi-square-1 tuning or post-result critical-value selection.

For p in [0,1], g(x)=p is the union of four rays with directions
v=(s1*p^2,s2*(1-p)^2), s1,s2 in {-1,1}. The minimum squared Mahalanobis
distance to a ray is d'Md-max(v'Md,0)^2/(v'Mv), M=V^-1. Invert its minimum
over signs. The full-rank ellipse excluding the origin is path-connected;
g is continuous there, so the image closure is an interval. Bracket each
endpoint between g(d) and 0 or 1; axis intersection supplies exact endpoints.
An ellipse containing the origin yields [0,1]; both-zero population contrasts
leave the ratio undefined, not equal to .5. Reject invalid/non-SPD covariance.

If the joint region covers the true nonzero contrast vector, its projection
covers g(d0). Thus coverage is at least joint coverage asymptotically under
iid sampling, correct working models, nonsingular contrast covariance, and
consistent observed-row inference. This is not an exact finite-sample claim
or a proof of uniform validity near simultaneous zero contrasts.

The geometry is inspired by von Luxburg & Franz (2009), Statistica Sinica
19:1095–1117, Sections 1.1, 2.1–2.2 (pp1097–1103), read directly at
https://www3.stat.sinica.edu.tw/sstest/oldpdf/A19n312.pdf . Their exact normal
ratio construction uses a one-dimensional Student-t critical value; our
two-dimensional Wald projection is deliberately different and potentially
conservative. Their theorem does not establish our fitted-contrast coverage.

## Fresh confirmation and unchanged judgment

Use the same seven DGP cells/truth integrals as the frozen boundary-correction
study: Gaussian/binary n=500/1000 regular, binary n=500 zero-direct,
zero-indirect, weak-direct. Each has 2,000 datasets. Seed = 20461020 +
10000*cell + replicate, disjoint from 20261020 and 20361020 studies.

Gaussian controls retain existing inference. Binary rows compare raw/Wald,
raw/projected, and previously frozen k=4/projected. All retain the existing
observed-row delta SE as a separate diagnostic; interval width is not a new
SE. Projection does not repair point bias. Record joint-region coverage,
contrast truth/estimates/covariance, projected limits, width, full-range rate,
and all errors. Report both exact-boundary stress results, not only regular
cells. Both-zero truth is tested for undefinedness, not assigned coverage.

Frozen diagnostics: zero failures; absolute bias/empirical SD <.10;
mean SE/empirical SD in [.90,1.10]; interval coverage in [.92,.98]. Apply and
report without rounding or changing limits, including conservative failures.
Report bias and coverage Monte Carlo SEs, raw replicates, versions, elapsed
time, source commit and source hashes. Independent tests compare projection
against dense ellipse boundaries and axis/scale/sign/swap identities; check
contrast covariance calibration and reproduce summary statistics separately.

Passing this finite study alone does not prove general boundary validity.
Any failed frozen cell retains both PR holds; full commit-bound package CI
is still required separately before a future merge/status transition.
