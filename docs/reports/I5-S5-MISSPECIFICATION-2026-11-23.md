# I5: Gaussian working-model misspecification in s5 at N = 1,000

## Targets

The s5 generator uses \(X\sim\mathrm{Bernoulli}(0.5)\), \(M=0.5X+U\), and
\(Y=0.1X+BM+E\), where \(U,E\) are independent standard normals and the
independent responder slope \(B\) is 0.8 with probability 0.7 and zero
otherwise. For the scenario's potential outcomes, the same \(B\) applies
to both arms, while the mediator and outcome innovations are drawn
independently across arms. Conditional on \(B=b\), the arm difference is
normal with mean \(0.5b\) and variance \(2(1+b^2)\). The mixture target is
therefore \(0.7\Phi(0.4/\sqrt{3.28})+0.3/2=0.5611803265\).
Numerical integration over the mediator difference confirmed the result
to within \(10^{-12}\).

The population Gaussian regression has mediator slope \(a=0.5\),
outcome slope \(b=E[B]=0.56\), and mediator residual variance 1.
Its outcome residual is \((B-E[B])M+E\), with marginal variance
\(1+\operatorname{Var}(B)E[M^2]=1.1512\). Substituting the working-model
moments in the Gaussian superiority formula yields \(0.5649726135\),
or \(0.0037922870\) above the mixture target. This difference persists
with exact Gaussian evaluations and arbitrarily large samples.

## Design

Fresh, paired comparison using the checked-in s5 generator and paper
MBCO routine from source commit `16bc700`. There were 2,000 independent
datasets of size 1,000. Each dataset was used for the MBCO interval and
the percentile interval from 1,000 pairs resamples. The bootstrap
evaluated the Gaussian plug-in score analytically in every resample,
eliminating finite-pair Monte Carlo error in that score. Both procedures
used nominal 95% intervals. Each set of endpoints was scored against
both targets, with missing endpoints counted as misses. The completed
run used R 4.6.2 and L'Ecuyer-CMRG seed 20261123 with a distinct stream per
dataset. The historical 0.92 result is the comparison motivating I5.
This run uses fresh streams.

From the repository root:

```sh
Rscript docs/reports/audit_i5_s5.R /tmp/probmed-i5-replications.csv 2000
```

The [replication records](data/i5-s5-2026-11-23/replications.csv.gz)
have uncompressed SHA-256
`cb229eeaef18f98d138aa0aa0d744c6a9bfad35878168980a2279289961676a2`.
An independent CSV audit checked the 2,000 unique keys, finite ordered
endpoints, matching pilot rows, and coverage counts. The point estimates
from the two implementations agreed within 1e-12 in every dataset.

## Results

| Method | Mixture target coverage (95% MC CI) | Gaussian target coverage (95% MC CI) | Paired gain (MCSE) | Mean width |
|---|---:|---:|---:|---:|
| MBCO | 1,852/2,000 = 92.6% (91.4–93.7%) | 1,894/2,000 = 94.7% (93.6–95.6%) | 2.1 pp (0.51 pp) | 0.03451 |
| Exact Gaussian percentile bootstrap | 1,854/2,000 = 92.7% (91.5–93.8%) | 1,896/2,000 = 94.8% (93.7–95.7%) | 2.1 pp (0.51 pp) | 0.03485 |

The coverage intervals are exact binomial Monte Carlo intervals.
Both methods returned finite intervals for all 2,000 attempts.
The common estimator had mean bias +0.003811 against the mixture target
(MCSE 0.000200) and +0.000019 against the Gaussian target. The sampling SD
was 0.008944 and the mean bootstrap SE was 0.008925. The analytic target
gap is 0.424 sampling SDs. A symmetric normal approximation with that
shift predicts 92.92% coverage.

Switching the target recovered 74 MBCO hits and lost 32, and recovered
73 bootstrap hits and lost 31. The direction agrees with the positive
bias. The mixture target lay below the lower endpoint 130 times for MBCO
and 124 times for the bootstrap, versus 18 and 22 misses above the upper
endpoint. Against the Gaussian target the corresponding splits were
56/50 and 51/53.

The paired comparison identifies a 2.1 percentage point loss in each
method from targeting the mixture with Gaussian working-model intervals
in this design. The undercoverage persists when the bootstrap score is
evaluated exactly. These findings apply to s5 at N = 1,000.

## Proposed manuscript wording for author review

> In the heterogeneous s5 scenario at N = 1,000, the Gaussian working-model
> limit exceeded the mixture estimand by 0.0038. In a fresh paired audit,
> coverage of the mixture estimand was 92.6% for MBCO and 92.7% for the
> bootstrap with exact Gaussian score evaluation. The same intervals
> covered the Gaussian working-model target 94.7% and 94.8% of the time,
> respectively, consistent with a coverage loss from misspecification.
