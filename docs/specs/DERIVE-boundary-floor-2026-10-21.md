# G2: minimax floor and weak-signal ladder (2026-10-21)

Plan: `tasks/plan-boundary-point-risk.md` (Phase G). Script: `scripts/check-boundary-floor-ladder.R`;
output: `docs/specs/data/check-boundary-floor-ladder-2026-10-21.txt`. No candidate estimator and no
archived estimate is used. Nothing here is a registered gate; the author's tolerance (G1) is separate.

## 1. Floor (Le Cam two-point)

For any estimator T and two data laws with truths theta0, theta1 (gap D = |theta1 - theta0|) and
total variation tau between the n-sample laws:

    max_i E_i (T - theta_i)^2 >= D^2 (1 - tau) / 8,   so   max RMSE >= D sqrt((1 - tau) / 8).

Derivation: let psi pick the nearer truth. Each wrong pick costs at least (D/2)^2, so
R0 + R1 >= (D/2)^2 (P0(psi = 1) + P1(psi = 0)) >= (D/2)^2 (1 - tau), and max >= average.
This is a lower bound on risk, not an attainable rate. It applies to any estimator, randomized or not.

## 2. Ladder at n = 500 (exact likelihood-ratio TV, Monte Carlo R = 20,000 datasets per row)

Positions are set by data-law TV to the boundary law, not by SE units of d2 (audit finding).

| Path | Boundary truth | TV target | Parameter at target |
|------|---------------:|----------:|--------------------:|
| direct (alpha .7, interaction 0) | P = 1 at direct 0 | .05 / .10 / .20 / .40 | .0175 / .0352 / .0713 / .1477 |
| indirect (direct .6, interaction .2) | P = 0 at alpha 0 | .05 / .10 / .20 / .40 | .0080 / .0160 / .0324 / .0677 |

Positions come from linear interpolation of the parameter on TV between grid rows, so they are
approximate; the unrounded table is in the output file.

## 3. Floors

| Path | Row (param) | Gap D | TV | Floor on max RMSE |
|------|------------:|------:|---:|------------------:|
| direct | .01 | .117 | .029 | .041 |
| direct | .05 (cell 7) | .228 | .142 | .075 |
| direct | .15 | .339 | .406 | **.092** (largest) |
| direct | .25 | .399 | .625 | .086 |
| indirect | .10 | .285 | .571 | **.066** (largest) |

**Reading.** Over the direct ladder no estimator can have worst-case RMSE below about .09 at
n = 500; over the indirect ladder, about .07. A per-cell gate tighter than that is infeasible by
construction, which is the same kind of finding as the bias/SD infeasibility in the T3 SPEC. The
floor is largest at middle TV (D grows while 1 - tau shrinks), not at the cusp.

## 4. Controls (script stops if any fails; all passed)

- shift 0 gives TV 0 on both paths.
- direct .05 gives TV .1418 vs the independent .1404 (`check-gate-feasibility-tv.R`), within 4 MC SE.
- largest shifts give TV > .9 on both paths (.9997 and 1.0000).

## 5. Limits and what this does not settle

- Boundary-to-point pairs only. Ladder-to-ladder pairs, and the max over the ladder of a combined
  risk, need the law at each point; not computed, so the floors are a lower envelope.
- The cusp is steep: D is .117 at direct .01 and .151 at alpha .02, so small shifts already carry
  large truth gaps. A tolerance must be read against D, not against the shift.
- Floors are for n = 500. Rate in n is an open assumption: Q5 was reversed, so n = 1000 needs its
  own floor from the same script with n changed, and its own tolerance.
- The tolerance and multiple (Q3) remain the author's G1 input; this document only says which
  tolerances are infeasible.

## 6. n = 1000 (added 2026-10-21; `Rscript scripts/check-boundary-floor-ladder.R 1000`)

Output: `docs/specs/data/check-boundary-floor-ladder-n1000-2026-10-21.txt`. Same laws, same controls
(the .1404 control applies at n = 500 only; at n = 1000 the cell 5 to 7 TV must exceed the n = 500 value, and does: .1981).

| Path | n = 500 largest floor | n = 1000 largest floor | Ratio |
|------|----------------------:|-----------------------:|------:|
| direct | .092 (param .15) | .082 (param .10) | .89 |
| indirect | .066 (param .10) | .059 (param .05) | .90 |

Ladder positions at TV .05 / .10 / .20 / .40:
direct .0125 / .0250 / .0505 / .1049; indirect .0057 / .0114 / .0231 / .0474
(each about .71 times its n = 500 position, consistent with 1/sqrt(n) scaling of the parameter at fixed TV).

**Reading.** Doubling n lowers the worst-case two-point floor only by about 10% (.092 to .082, .066 to .059),
not by the factor .84 that a pure n^(-1/4) rule would give or .71 for n^(-1/2). The floor is a maximum over a
ladder whose TV-matched positions move with n, so it does not scale like either rate; a tolerance at n = 1000
must come from its own floor, as the reversal of Q5 requires, and cannot be derived from the n = 500 tolerance.

## 7. In-between and node-to-node pairs (G2b, added 2026-10-21)

Script `scripts/check-boundary-floor-pairs.R [n]`, outputs `docs/specs/data/check-boundary-floor-pairs-n{500,1000}-2026-10-21.txt` and
`floor-pairs-n{500,1000}.csv` (all pairs). Grid per path: boundary, the four TV-matched positions, the cell-7-like .05, and midpoints
of consecutive nodes (R = 10,000 datasets per base law). Each pair uses data from law i and the exact likelihood ratio of law j.
Controls passed: identical-law TV 0, and the boundary pair at .05 reproduces the section 2 TV within 4 combined MC SE at both n.

| n | Path | Largest boundary-to-node floor | Largest node-to-node floor |
|---|------|-------------------------------:|---------------------------:|
| 500 | direct | .0926 | .0634 (.0088 vs .1477) |
| 500 | indirect | .0675 | .0481 (.0040 vs .0677) |
| 1000 | direct | .0820 | .0568 (.0063 vs .1049) |
| 1000 | indirect | .0591 | .0429 (.0029 vs .0500) |

**Reading.** The boundary pair binds: node-to-node floors are about 70% of it, because the truth gap between two interior points
is smaller while the TV between them is also smaller. So the floors in sections 3 and 6 stand as the two-point floor, and an
in-between control (G1c) does not raise the floor. It still matters as a guard against estimators that pass only at the nodes.

**Limits.** These are two-point bounds. A multi-point (Fano-type) bound over the whole ladder can be larger; it is not computed here,
so the floor remains a lower bound on max RMSE, not a tight one. MC SE is in the output; the boundary maxima differ from sections 3
and 6 by up to .0003 because R here is 10,000 and the seed differs.
