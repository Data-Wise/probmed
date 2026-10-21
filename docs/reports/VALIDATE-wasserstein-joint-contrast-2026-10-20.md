# Joint-contrast confidence-set investigation

**Decision: keep PRs #54/#55 blocked and in draft.** Projection substantially
improves boundary coverage but does not pass the frozen diagnostics. No public
estimator, inference fallback, required status, or hosted workflow changed.
Full exact-commit software CI remains pending; it cannot override these failures.

## Frozen method and evidence scope

The [plan](PLAN-wasserstein-joint-contrast-2026-10-20.md) was committed as
549daee before implementation. Source f5669a6 was committed before generating
fresh results. The research prototype projects a 95% two-dimensional Wald
ellipse for the signed binary probability contrasts through the square-root
absolute-ratio map. It includes both contrast signs and their covariance.
The origin is excluded from the ratio's domain; an ellipse containing it returns
the full [0,1] range, not an estimate of .5. Both population contrasts zero is
undefined and is not a coverage cell. This method is binary-specific; Gaussian
controls retain the previously implemented inference.

For a nonzero true contrast vector, joint coverage implies projected-ratio
coverage by set inclusion. Under an appropriate joint CLT and consistent
nonsingular covariance, this yields an asymptotically conservative procedure
without differentiating the ratio at a zero component. It does not establish
finite-sample exactness, uniform validity near simultaneous zero, correct
working-model specification, or validity of the public bootstrap fallback.

Primary source: von Luxburg & Franz (2009), *Statistica Sinica* 19:1095–1117,
[Sections 1.1 and 2.1–2.2, pp1097–1103](https://www3.stat.sinica.edu.tw/sstest/oldpdf/A19n312.pdf),
read directly. Their geometric treatment explains why ratio confidence sets
must allow weak-denominator uncertainty. Their exact normal-ratio construction
uses a one-dimensional Student-t critical value, unlike our joint chi-square-2
Wald region. Their theorem is motivation, not a validation theorem for our
fitted contrasts or nonlinear absolute-ratio mapping.

## Fresh simulations: all frozen cells retained

Seven cells, 2,000 datasets each, seed base 20461020; the DGP/truth integrals
are unchanged from the earlier frozen boundary study. All 14,000 datasets
finished with zero failed fits, yielding 34,000 paired method rows in 274.20
seconds under R 4.6.1. Original raw estimates, raw/projected, and previously
frozen k=4/projected were compared for binary outcomes. The corrected point
estimate does not change the projected interval. Existing delta SEs remain
separate diagnostics, not SEs inferred from projected width.

| Cell | Point/interval | Bias/SD | SE/SD | Coverage | Frozen result |
| --- | --- | ---: | ---: | ---: | --- |
| Gaussian n=500 | unchanged | .0413 | .9894 | .9485 | Pass |
| Binary n=500 | k=4/projected | .0696 | 1.0260 | .9800 | Pass |
| Gaussian n=1,000 | unchanged | .0592 | 1.0156 | .9535 | Pass |
| Binary n=1,000 | k=4/projected | .0359 | 1.0057 | .9885 | Fail: coverage |
| Binary zero direct | k=4/projected | 3.2153 | 3.0630 | .9875 | Fail: all three |
| Binary zero indirect | k=4/projected | 2.6022 | 2.7112 | .9910 | Fail: all three |
| Binary weak direct | k=4/projected | 1.0577 | 2.7541 | .9890 | Fail: all three |

Limits remain exactly bias/SD <.10, SE/SD [.90,1.10], coverage [.92,.98],
and zero failures. Coverage .9800 passes the inclusive upper limit; .9885
does not. Raw binary n=500 still fails bias/SD (.1236); its coverage is .9560
with Wald and .9800 with projection. No candidate was selected or retuned
using these results.

| Binary cell | Raw Wald coverage | Projected coverage | Coverage MCSE | Mean projected width | Width MCSE | Joint-vector coverage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Regular n=500 | .9560 | .9800 | .00313 | .4252 | .00306 | .9415 |
| Regular n=1,000 | .9610 | .9885 | .00238 | .2579 | .00216 | .9410 |
| Zero direct | .4875 | .9875 | .00248 | .5512 | .00116 | .9465 |
| Zero indirect | .5510 | .9910 | .00211 | .6269 | .00574 | .9445 |
| Weak direct | .8025 | .9890 | .00233 | .5537 | .00124 | .9435 |

The zero-indirect interval is the full [0,1] in 29.65% of datasets; no other
cell has a full-range interval. All raw and k=4 point estimates lie inside
their reported intervals in this study. Mean estimated contrast covariance
relative to empirical covariance has generalized eigenvalues .9483–1.0865
across the binary cells. Joint-vector coverage is .9410–.9465, not a claim
of exact .95 calibration. Its binomial MCSE is approximately .0050–.0053.
Projection is intentionally more conservative than joint-vector coverage.

## Why this does not unblock merging

In a regular scalar direction, projecting the two-dimensional ellipse gives
first-order scalar coverage 2*Phi(sqrt(qchisq(.95,2)))-1 = .9856247, not .95.
At a true zero component with the other fixed nonzero, covering that endpoint
also corresponds to the ellipse intersecting the relevant axis, yielding the
same limiting marginal-normal probability. Thus coverage above the frozen
.98 ceiling is structurally plausible, not something to round away. These
approximations presume regular joint-contrast inference and nonsingularity.

More importantly, a new interval does not repair point-estimate bias or turn
the nonregular ratio into a regular delta-method estimator. The bounded k=4
correction still fails bias and SE diagnostics in every stress cell. We do
not reinterpret those failures as passing because interval coverage improved.
The numerical study supports a useful conservative confidence-set prototype,
not a merge-ready general boundary repair.

Next research should investigate boundary point estimation and nonregular
sampling uncertainty jointly, with independently justified calibration for
test inversion if pursued. Any modified estimator or critical value requires
a new frozen design and independent fresh seeds, not tuning this set.
No frozen diagnostic or scientific hold is relaxed here.

## Reproducibility and verification

- `Rscript scripts/test-wasserstein-joint-contrast.R`: 19 assertions pass.
  An independent 200,001-angle ellipse oracle covers signed/axis-crossing
  fixtures; scale, swap, sign and invalid-covariance tests also pass. Replacing
  the procedure in memory by [0,1] fails all four dense-boundary comparisons.
- Another 1,000 randomly generated inside-ellipse points (seed 814) satisfy
  mapped containment. This is numerical evidence, not a symbolic proof.
- `Rscript scripts/validate-wasserstein-joint-contrast.R OUTPUT_DIRECTORY`
  reproduces the prespecified simulation. Each binary cell's first ten fits
  independently agree with production points/SEs to 1e-9/1e-8.
- `Rscript scripts/audit-wasserstein-joint-contrast.R` independently recomputes
  all 17 saved summary rows and verifies 14,000 distinct seeds with no overlap
  with either earlier study, all error fields empty, joint quadratic-form
  coverage, and joint-to-projected containment for every saved binary row.
- [Raw results](data/wasserstein-joint-2026-10-20/replicates.csv),
  [summary](data/wasserstein-joint-2026-10-20/summary.csv), and
  [provenance](data/wasserstein-joint-2026-10-20/provenance.txt) retain all
  estimates, contrasts, covariance entries, intervals, truths, seeds and MCSEs.
  Replicates SHA-256:
  `fc6a4d8745f448890642a1ef169c02322fd71e4697d25faaf73df4778fd6bc34`.
  Projection source SHA-256:
  `fe64fe60b3db0fa81e4475ef529211460bd2c6579a8507703711eb68ddc44afe`.

The prototype shares the earlier binary IF implementation and quadrature nodes;
its geometry oracle is independent, but this is not an independent refitting
implementation or proof of all model assumptions. No package API changed,
so this research-only check is not a substitute for full package CI or
public-path inference validation. Existing I6/I7/release/submission holds persist.
