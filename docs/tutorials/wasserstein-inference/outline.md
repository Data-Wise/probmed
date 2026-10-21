# When good coverage is not enough

## Proposed technical mini-book: point estimation and confidence sets near mediation boundaries

Status: **proposed architecture only**, not written chapters or a rendered book.
Requested stage: Tome outline. Native preflight is blocked: PyYAML is absent
from the checked Python runtimes and this package has no `.flow/research-config.yml`.
No package-root research configuration, notebook, vault, or rendering engine
has been installed or silently created. This outline is a review artifact,
not a claim that the configured Tome planning gate passed.

### Thesis, audience, and teaching contract

Thesis: uncertainty must follow the observed data that were fitted, and a valid
confidence set does not automatically supply a low-bias point estimator or a
calibrated standard error at a nonsmooth boundary.

Assumed floor: college algebra, introductory differentiation/integration, and
basic probability. Teach matrix covariance, asymptotic notation, estimating
equations, quantiles, test inversion, and simulation diagnostics before using
them. An optional graduate route skips primers, not assumptions or caveats.

Proposed `book.register`: **tutorial**. Proposed `voice.register`: **methods**.
House style: US English; explicit assumptions; plain-language meaning before
symbols; no unexplained algebra; distinguish a derivation, a computational
check, and empirical evidence. These are proposals, not persisted configuration.

Recommended spine: **A, outline-first**. Verified authored-source inventory:
one repository, zero Obsidian sources, zero separate teaching sources. The
argument is recoverable from existing code, frozen plans, and archived studies.
External papers are supplemental research, not additional authored repositories.
No existing source is superseded; reports and source hashes remain authoritative.

Every chapter will contain objectives, prerequisite primers, a geometry/intuition
view where honest, full derivations, at least two worked examples for substantial
results, faded examples, self-explanation prompts, checkpoints with solutions,
misconception checks, and an end-of-chapter exercise set with a separate key.
Nothing below claims those devices have already been authored or verified.

## Authored-source map

Paths below resolve from the probmed worktree root, three directories above
this directory. Stable source identifiers are used throughout the architecture.

| ID | Source | Contribution |
|---|---|---|
| S01 | `AGENTS.md`; `R/wasserstein-pmed.R` | Estimand distinction, corner laws, public API versus observed-corner engine |
| S02 | `R/wasserstein-observed.R` | Fitted-model g-computation, quadrature, empirical-covariate terms and row influence functions |
| S03 | `docs/reports/REVIEW-dev-main-blocker-2026-10-20.md` | Simulated-draw SE bug and fixed-data counterexample |
| S04 | `docs/reports/VALIDATE-wasserstein-observed-inference-2026-10-20.md`; `scripts/validate-wasserstein-observed.R` | Original four-cell fresh validation and frozen diagnostics |
| S05 | `docs/reports/INVESTIGATE-wasserstein-binary-bias-2026-10-20.md`; `scripts/diagnose-wasserstein-binary-bias.R` | Exact descriptive bias decomposition and replay |
| S06 | `scripts/pilot-wasserstein-bias-correction.R`; `scripts/diagnose-wasserstein-jackknife-extreme.R` | Grouped/delete-one jackknife failure evidence |
| S07 | `scripts/wasserstein-boundary-candidates.R`; `scripts/test-wasserstein-boundary-candidates.R` | Curvature correction, independent contrast implementation, derivative checks |
| S08 | `docs/reports/VALIDATE-wasserstein-boundary-correction-2026-10-20.md` | k=4/k=9 regular success and boundary failure |
| S09 | `scripts/wasserstein-joint-contrast.R`; `scripts/test-wasserstein-joint-contrast.R` | Ellipse-to-ratio projection and geometry oracle |
| S10 | `docs/reports/VALIDATE-wasserstein-joint-contrast-2026-10-20.md`; `scripts/audit-wasserstein-joint-contrast.R` | Joint versus scalar coverage, bias/SE failures and independent audit |
| S11 | `scripts/wasserstein-boundary-selection.R`; `scripts/test-wasserstein-boundary-selection.R` | Hard selection, constrained center, Gaussian proxy SE, Fieller inversion |
| S12 | `docs/reports/PLAN-wasserstein-boundary-selection-2026-10-20.md`; `docs/reports/VALIDATE-wasserstein-boundary-selection-2026-10-20.md` | Latest frozen design, numerical amendment, fresh findings and limitations |
| S13 | `scripts/validate-wasserstein-boundary-selection.R`; `scripts/audit-wasserstein-boundary-selection.R` | Driver reuse, paired comparisons, seed audit and direct linear-null oracle |
| S14 | `docs/reports/data/wasserstein-{observed,boundary,joint,selection}-2026-10-20/` | Four separately frozen raw evidence archives; braces denote four actual directories |
| S15 | `.STATUS` | Scientific holds and distinction from exact-commit software CI |

Supplemental research already read for the investigation:

- von Luxburg and Franz (2009), *A Geometric Approach to Confidence Sets for
  Ratios*, Statistica Sinica 19, 1095–1117; actual pages 1097–1105,
  sections 1.1, 2.1–2.2 and 3.2. Source:
  https://www3.stat.sinica.edu.tw/sstest/oldpdf/A19n312.pdf.
  Supports ratio-test geometry, not exact inference for our fitted contrasts.
- Pötscher and Schneider (2010), arXiv:0806.1652v3; actual introduction and
  sections 2–3 through Proposition 2, pages 1–4. Source:
  https://arxiv.org/pdf/0806.1652. Supports thresholding/nonuniformity cautions,
  not validation of this ratio or its Gaussian proxy SE.

Reference gap: DOI/Zotero-entry and attachment-key inventories are not yet
settled for a book bibliography. Do not mark metadata-only leads as read.

## Dependency-ordered chapter architecture

Stable labels below own chapter identity; printed numbers are positions only.

### 1. Which quantity are we estimating? `\label{ch:estimands}`

**Objective:** distinguish causal meaning from the numerical functional being
estimated. **Prerequisites:** assumed floor. **Sources:** S01, S15.

- Define observed treatment A, mediator M, outcome Y and baseline C; introduce
  potential outcomes and the corner law nu(a,a')=Law(Y(a,M(a'))).
- State consistency, positivity, relevant confounding/cross-world assumptions,
  and working-model scope separately; marginal laws are not a complete coupling.
- Distinguish the package's tie-adjusted stochastic-superiority P_med from the
  separate transport score P_med^W = W_I/(W_I+W_D). These are not synonyms.
- Define W_I=W2(nu11,nu10), W_D=W2(nu10,nu00), W_T=W2(nu11,nu00).
  Explain the triangle gap W_I+W_D-W_T, not an additive causal decomposition.
- Plan paired examples where transport distances discard direction but a
  superiority probability retains it. Explain why binary support is not evidence
  of ordinal/count support and a null score is not absence of every mechanism.
- Exercise: identify claims that require a joint counterfactual coupling versus
  marginal corner laws. Answer key must preserve these distinctions.

### 2. Transport geometry from the ground up `\label{ch:transport}`

**Objective:** derive the functional actually computed. **Depends on:**
`ch:estimands`. **Sources:** S01, S02.

- Primer: CDFs, quantiles, coupling, and squared transport cost.
- Derive W2^2(F,G)=integral_0^1 [Q_F(u)-Q_G(u)]^2 du in one dimension,
  stating finite-second-moment and monotone-coupling conditions.
- Work through normal laws: W2^2=(mu1-mu2)^2+(sigma1-sigma2)^2.
- Derive Bernoulli W2^2=|p-q| by matching common mass and moving unmatched
  probability over distance one. Do not substitute |p-q| for W2 itself.
- Let d1=p11-p10, d2=p10-p00. Obtain
  g(d)=sqrt(|d1|)/(sqrt(|d1|)+sqrt(|d2|)); map axes to P=0 and P=1.
- Worked fixtures include sign reversal, equal components, a single zero, and
  both zeros. Prove that double zero is undefined, not automatically 1/2.

### 3. Two sample sizes, two kinds of uncertainty `\label{ch:uncertainty}`

**Objective:** understand the original scientific blocker. **Depends on:**
`ch:transport`. **Sources:** S01–S04.

- Separate observed sample size n, inner simulation count K (public n_mc),
  bootstrap resamples B_boot, and validation repetitions R. Never call all B.
- Primer: conditional versus repeated-sampling variance and total variance.
  Specify what is held fixed in each experiment.
- Explain the legitimate directly observed-corner experiment and why its
  potential-based IF does not automatically cover corners simulated from fits.
- Reproduce the fixed-500-row counterexample: increasing draws 2000 to 200000
  changes the legacy SE .03781 to .003726 without creating new observations.
- Derive why Monte Carlo noise can vanish while fitted-model uncertainty remains.
- Checkpoint: distinguish an invalid uncertainty calculation from an inaccurate
  numerical approximation. Explicitly separate the legacy engine from public repair.

### 4. Deterministic fitted-model g-computation `\label{ch:gcomp}`

**Objective:** reconstruct every corner from observed data. **Depends on:**
`ch:uncertainty`. **Sources:** S02, S07.

- Write the Gaussian mediator regression and Gaussian/logistic outcome regression,
  including A-by-M interaction, model matrices, residual variances and C averaging.
- Derive binary corner means as Gaussian expectations of expit(eta); explain why
  expit(E eta) is not E expit(eta).
- Teach Gaussian-Hermite nodes/weights and the eigenvalue construction; distinguish
  integration accuracy from statistical sampling accuracy.
- Derive Gaussian conditional corner mixtures, numerical CDF inversion and
  Gauss-Legendre integration over quantile probability u.
- Compare quadrature with adaptive integration as an independent numerical oracle.
- Work through family/link scales, model-rank checks, convergence, positive variance,
  weights, root bracketing and deterministic reproducibility.

### 5. Influence functions and the sandwich `\label{ch:influence}`

**Objective:** propagate all fitted-observation uncertainty. **Depends on:**
`ch:gcomp`. **Sources:** S02, tests in `tests/testthat/test-wasserstein-observed.R`.

- Primers: vectors/matrices, covariance, Taylor expansion, estimating equations,
  asymptotic linearity, and the central limit theorem.
- Derive parameter-row IFs: Q_M^{-1} X_M residual_M; logistic
  Q_Y^{-1} X_Y(Y-p); and residual^2 minus estimated residual variance under
  the stated correctly specified working-model setting.
- Stack row contributions before computing covariance; the same observed rows
  generate corner estimates, so independent-corner variance addition is wrong.
- Derive binary corner IF = gradient_theta p times parameter IF plus the empirical
  covariate averaging term h(C_i)-mean h(C). Explain the omitted-C mutation.
- Derive Gaussian quantile IF by implicit differentiation of F(q)=u, including
  its density denominator and empirical CDF contribution; then transport IF.
- Chain to IF_P=[(1-P)IF_WI-P IF_WD]/(WI+WD), and SE=sd(IF_P)/sqrt(n).
- Explain why a working-model row IF is not automatically a semiparametric
  efficient influence function; the public string `eif` is not a proof.

### 6. Validation is a set of questions `\label{ch:diagnostics}`

**Objective:** interpret the frozen contract without substituting one metric
for another. **Depends on:** `ch:influence`. **Sources:** S04, S12–S14.

- Define bias, empirical SD, estimated SE, coverage, width, RMSE and fit failure.
- Derive MCSE(bias)=empirical SD/sqrt(R) and MCSE(coverage)=sqrt(c(1-c)/R).
  Discuss why uncertainty in ratios and tail diagnostics is more involved.
- Retain gates: zero failures; |bias|/SD<.10; mean SE/SD in [.90,1.10];
  coverage in [.92,.98]. These are frozen engineering diagnostics, not universal
  statistical laws or confidence bounds proving population bias below .10 SD.
- Work through binary n=500: bias/SD=.142817, SE/SD=.988387, coverage=.9630.
  Explain why this still fails. Separate inconclusive precision from gate failure.
- Distinguish these Wasserstein diagnostics from the separate exact-Gaussian
  P_med coverage gate, I6, and software CI. No cross-task gate substitution.
- Exercise: construct a wide interval with good coverage and a poor point estimate.

### 7. Locating the bias before correcting it `\label{ch:bias-replay}`

**Objective:** separate upstream fit bias from nonlinear transformation behavior.
**Depends on:** `ch:diagnostics`. **Sources:** S05, S14.

- Derive the exact finite-replicate identity:
  mean g(p_hat)-g(p_true) = [g(mean p_hat)-g(p_true)] +
  [mean g(p_hat)-g(mean p_hat)]. Define p as the three-corner vector here.
- Reproduce 2000 saved estimates within 1e-8; audit sign reversals and row identity.
- Explain .009851444 total bias, .001396391 upstream component and .008455053
  transformation component (85.8%). A decomposition is not a causal attribution
  theorem, a proof of bias order, or an immediately usable correction.
- Numeric exercise: the same identity for a nonlinear scalar function with
  known inputs; distinguish observed replay identities from population expectations.

### 8. Why jackknife extrapolation failed `\label{ch:jackknife}`

**Objective:** derive the method and diagnose its instability. **Depends on:**
`ch:bias-replay`. **Sources:** S06 and the completed pilot section of S05.

- Derive bias cancellation from E T_n=theta+a/n+smaller terms, first for
  delete-one and then for 20 equal delete-groups: 20 T_full-19 mean T_minus_group.
- State the smooth-expansion assumption; grouping is not an interchangeable
  substitute for an independently validated delete-one procedure.
- Reproduce the 200-development-dataset pilot: RMSE .067798 to .111977;
  group-corrected selected estimate .868992 to 1.913700; ordinary delete-one
  on that dataset gives -1.764091.
- Teach cancellation, leverage, instability near a cusp, and amplification by
  extrapolation. Explain why clipping changes the estimator and its bias.
- Explain why reusing an old SE for a corrected estimator is only exploratory.
- Exercise: find the invalid step in applying smooth bias cancellation at a cusp.

### 9. Bounded curvature correction `\label{ch:curvature}`

**Objective:** derive the k=4/k=9 research candidates completely. **Depends on:**
`ch:bias-replay`, `ch:influence`. **Sources:** S07–S08.

- Primer: multivariate Hessian, trace, delta-method curvature, logit and expit.
- At nonzero signed d1,d2, derive g1=P(1-P)/(2d1),
  g2=-P(1-P)/(2d2), and all three Hessian entries; verify by independent finite
  differences rather than claiming a symbolic backend was used.
- Derive second-order bias .5 trace(H_g V). Translate its first-order correction
  into logit shift; distinguish correcting transformation bias from all bias.
- Derive the implemented regularization: D_j=d_j^2+k V_jj,
  r_j=V_jj/D_j, r12=V12 d1 d2/(D1 D2),
  s=-(1+2P)r1/8+(3-2P)r2/8-(1-2P)r12/4,
  corrected P=expit(logit(P)-s).
- Explain endpoint handling, undefined double zero, scale/sign/swap identities,
  and how a PSD covariance controls the damped shift.
- Reproduce fresh regular success for k=4; explain k=9 bias/SD=.100275 is a
  failure, not .10 rounded into a pass. Boundary coverage .502/.5605/.7855
  exposes the limit of boundedness as a remedy.

### 10. Why boundaries change the mathematics `\label{ch:nonregular}`

**Objective:** understand why ordinary delta-method reasoning breaks.
**Depends on:** `ch:curvature`, `ch:transport`. **Sources:** S02, S08, S12;
supplemental thresholding paper for cautions, not a theorem about this score.

- Define regularity as used here, nondifferentiability, local alternatives,
  pointwise versus uniform statements, and convergence-rate notation.
- Derive, under a nondegenerate root-n contrast limit, why sqrt(abs(d_hat1))
  is order n^{-1/4} at d1=0 while a separated d2 stays fixed.
- Teach the resulting folded/square-root-normal limit and why a normal SE
  approximation is not automatically restored by a bounded transformation.
- Analyze paths where both contrasts approach zero; distinguish failure of
  concentration of a ratio from a ratio defined at population double zero.
- Derive the elementary bounded-estimator fact: T>=0 and E T=0 imply T=0
  almost surely. Do NOT turn this into an unproved impossibility claim for
  the actual nonzero bias-tolerance gates.
- Exercise: contrast fixed separated, exact-axis, and shrinking-signal sequences.

### 11. Projecting joint confidence regions `\label{ch:joint-region}`

**Objective:** learn a valid coverage argument and its conservatism.
**Depends on:** `ch:nonregular`, `ch:influence`. **Sources:** S09–S10.

- Primer: positive-definite covariance, Mahalanobis distance and chi-square df.
- Construct E={x:(x-d_hat)' V_hat^{-1}(x-d_hat)<=chi2_.95,2} and project its
  image closure through g. Prove coverage inheritance when the true contrast
  vector lies in E; label the joint CLT approximation correctly.
- Explain origin containment returning [0,1], axis intersections, four signed
  rays, interval closure, and numerical root bracketing.
- Derive fixed-p directions (p^2, plus/minus (1-p)^2) and the projection metric;
  compare against the independent dense ellipse-boundary oracle.
- Explain asymptotic scalar projection coverage
  2 Phi(sqrt(chi2_.95,2))-1=.9856247 in a regular scalar setting.
- Reproduce excessive coverage and persistent bias/SE failures: conservative
  interval projection is not a complete repair and can fail an upper coverage gate.

### 12. Fieller-style union-null test inversion `\label{ch:fieller}`

**Objective:** calibrate the hypothesis about the scalar functional directly.
**Depends on:** `ch:joint-region`. **Sources:** S11–S13; von Luxburg–Franz.

- Primer: confidence sets as inversion of tests; union versus intersection nulls.
- For p in [0,1], set a_s(p)=(s(1-p)^2,-p^2), s in {-1,+1}.
  Derive T_s(p)=(a_s' d_hat)^2/(a_s' V_hat a_s).
- Derive the equivalence between standardized perpendicular contrast and
  Mahalanobis distance to a full line. Show why sign lines collapse on axes.
- Invert min_s T_s(p)<=chi2_.95,1. Under a regular joint CLT and a true line,
  its rejection event is contained in rejection of that line's test; give the
  scoped asymptotic coverage argument, not exact finite-sample normality.
- Explain why this is not relabeling a 95% joint ellipse. Verify endpoint
  equations with an independent linear-null oracle.
- Reproduce fresh binary coverages .9475/.9500/.9440/.9450/.9715 and their MCSEs.
  Discuss possible conservatism near intersecting sign lines and undefined origin.

### 13. Hard boundary selection and constrained centers `\label{ch:selection}`

**Objective:** derive both point-selection rules and their failure mechanisms.
**Depends on:** `ch:fieller`, `ch:nonregular`. **Sources:** S11–S13.

- Define weak_j iff |d_hat_j|/sqrt(V_jj)<=lambda. Explain thresholds 1.959964
  and 2.807034 as fixed candidates, not posterior probabilities of a zero effect.
- Exactly one weak: choose its axis (P=0 for first, P=1 for second).
  Neither or both weak: retain k4. Record ambiguous both-weak cases.
- Derive the constrained center mu*=d_hat-V[,j] d_hat_j/V_jj by minimizing
  Mahalanobis distance subject to mu_j=0. Explain the correlation adjustment.
- Plan examples with correlated versus independent contrasts and changing axes.
- Show why selecting a boundary more often can reduce exact-null bias yet
  introduce false-axis bias under regular or weak alternatives.
- Fresh evidence: regular n=500 selection 22.95%/52.70%; zero-indirect ambiguous
  cases 15.45%/41.95%; weak-direct stronger-selection bias/SD=7.8773.
- Discuss point/interval compatibility: stronger-selector point-in-interval
  rates .7025–.9670; never conceal this by a post hoc interval hull.

### 14. Estimating uncertainty after selection `\label{ch:proxy-se}`

**Objective:** understand what the Gaussian-experiment proxy computes and omits.
**Depends on:** `ch:selection`. **Sources:** S11–S13; thresholding-paper cautions.

- Generate X*=mu*+Z L, with L' L=V for row-vector normal draws. Apply the
  entire selection/k4 rule, including switching, and take its Monte Carlo SD.
- Teach Cholesky orientation, antithetic pairs, fixed common random numbers,
  and the difference between simulation precision and statistical validity.
- Explain the initial 32768-draw precision failure (ratio 1.102004), the
  documented pre-dataset increase to 65536, and independent reference checks.
- Distinguish this constrained Gaussian proxy from an ordinary nonparametric
  bootstrap, a fitted-data parametric bootstrap, and a proof of uniform validity.
- Reproduce zero-direct stronger-selection bias/SD=.0761 but SE/SD=1.9622;
  primary regular n=1000 SE/SD=1.2818. Candidate selection changes uncertainty.
- Exercise: diagnose a rule that returns almost always the same wrong endpoint.

### 15. Building an auditable simulation study `\label{ch:simulation}`

**Objective:** reproduce all design, execution and evidence safeguards.
**Depends on:** `ch:diagnostics`, `ch:proxy-se`. **Sources:** S04, S08, S10, S12–S14.

- ADEMP: aims, data-generating mechanisms, estimands, methods and performance.
- Specify C~N(0,1), A|C~Bern(expit(.3C)),
  M=alpha A+.5C+epsilon_M, eta=-.5+direct A+.8M+interaction AM+.3C,
  with Gaussian or Bernoulli Y. Explain all seven cells and n=500/1000.
- Derive population corners by integrating the joint Gaussian linear predictor;
  derive independent adaptive-integration truth, including zero/weak scenarios.
- Teach paired method comparisons, independent datasets, disjoint seed bases,
  worker allocation, deterministic integrals, and failure retention.
- Distinguish development replay, targeted extreme diagnostics, and fresh
  confirmation. Never tune a candidate on its claimed independent validation set.
- Walk through 14000 datasets versus 44000 method rows, cell checkpoints,
  full CSV archives, source commits, hashes, runtime and R-session provenance.
- Audit CSV type conversion (`error` must stay character), all summary formulas,
  seed disjointness, selection centers and test-inversion inclusion independently.

### 16. Scientific conclusions and the next admissible question `\label{ch:synthesis}`

**Objective:** distinguish accomplishments from unresolved inference.
**Depends on:** every preceding chapter. **Sources:** S03–S15.

- Build an evidence ladder: derivation, numerical agreement, mutation detection,
  repeated-sampling validation, independent audit, public-path tests, full CI.
- Compare approaches without collapsing their different target problems:
  observed-row repair fixes omitted uncertainty; curvature correction helps
  regular bias; Fieller inversion helps interval calibration; hard selection
  does not jointly satisfy point-bias and SE diagnostics.
- Mark what was not established: general boundary-bootstrap coverage, uniform
  local-alternative validity, double-zero estimability, model-misspecification
  robustness, or a production-ready point correction.
- Frame next work as a risk/bias trade-off and smooth constrained-estimation
  problem. Potential candidates are research questions, not recommendations
  already validated. Preserve all existing cells, gates, and fresh-seed rules.
- Capstone: audit a deliberately flawed report claiming that acceptable coverage
  alone authorizes merging. Separate scientific hold, software CI and author authority.

## Appendices and reusable teaching assets

- `app:notation`: symbol ledger and bridge to literature; each symbol has one
  meaning per context. Distinguish regression coefficient beta from validation
  repetitions R, covariance V from residual variance, and p from P_med^W.
- `app:matrix-primer`: covariance, PSD/SPD, inversion, trace, Cholesky and
  quadratic forms; each operation has a numerical check.
- `app:derivative-ledger`: complete gradients/Hessians, implicit quantile
  derivatives, influence-function chain rules and legal-step justifications.
- `app:proof-ledger`: claims with assumptions, proof dependencies, theorem versus
  conjecture labels, and counterexamples; no theorem inferred from a green test.
- `app:reproduction`: exact commands, immutable source revisions, seed maps,
  hashes and read-only scripts for each archived study.
- `app:solutions`: worked solutions, hint ladders, error-analysis exercises and
  cumulative final with a section-addressed self-diagnostic.

## Reading paths

- Full derivation path: chapters in the dependency order above.
- Inference-focused graduate path: `ch:estimands`, `ch:uncertainty`,
  `ch:influence`, `ch:diagnostics`, `ch:nonregular`, `ch:joint-region`,
  `ch:fieller`, `ch:selection`, `ch:proxy-se`, `ch:synthesis`.
- Forensic reproducibility path: `ch:diagnostics`, `ch:bias-replay`,
  `ch:jackknife`, `ch:curvature`, `ch:simulation`, then `ch:synthesis`.

## Definition ownership and pre-draft checks

No chapter supersedes an existing report. Future prose must introduce definitions
before uses, with stable labels owned as follows:

| Terms requiring definitions | Owner |
|---|---|
| Corner law, coupling, stochastic superiority, transport score | ch:estimands |
| Quantile, W2, signed contrast, axis, undefined double zero | ch:transport |
| Conditional variance, Monte Carlo error, observed-data uncertainty | ch:uncertainty |
| G-computation, quadrature, working model, empirical-covariate distribution | ch:gcomp |
| Estimating equation, IF, sandwich, CLT, asymptotic linearity | ch:influence |
| Bias, RMSE, coverage, MCSE, frozen diagnostic | ch:diagnostics |
| Replay identity, transformation bias | ch:bias-replay |
| Jackknife, extrapolation, leading bias expansion | ch:jackknife |
| Hessian, curvature, trace, damping, logit | ch:curvature |
| Nonregularity, local alternatives, pointwise/uniform validity | ch:nonregular |
| Mahalanobis metric, joint confidence region, image closure | ch:joint-region |
| Test inversion, Fieller set, union null, degrees of freedom | ch:fieller |
| Hard threshold, constrained center, false-axis selection | ch:selection |
| Gaussian proxy SE, antithetic draw, selection switching | ch:proxy-se |
| ADEMP, paired comparison, provenance, development/confirmation split | ch:simulation |

This is a prospective ownership audit, not a passed terminology-closure audit
of unwritten chapters. Before chapter drafting: settle local book configuration
and runtime, preserve the proposed floor/notation, and resolve source-map gaps.
Before rendering: run the actual closure/citation/derivation verification gates,
then inspect real HTML and PDF artifacts. No rendering is claimed here.
