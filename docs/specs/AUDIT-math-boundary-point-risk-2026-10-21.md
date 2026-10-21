<!-- Independent math audit by a fresh-context subagent (savant proof-architect L1-L4 method), 2026-10-21. Model output, not verified by the author; its findings were checked and applied selectively (see the commit message and the corrected docs). Treat as data. -->

# Math audit: boundary point-risk worktree (recomputed in R; scratch scripts a.R, tv.R, emp.R, d1.R, d2.R in this folder)

| id | claim | verdict | evidence / recomputed | sev |
|---|---|---|---|---|
| 1a | target P=sqrt|d1|/(sqrt|d1|+sqrt|d2|); truths cell5=1, cell6=0, cell7=.7721 | CONFIRMED | rebuilt from DGP script: d=(.11230,0) P=1; d=(0,.11908) P=0; d=(.11262,.009811) P=.772109 | - |
| 1b | E|Z|^(1/2)=2^(1/4)G(3/4)/sqrt(pi)=.822 | CONFIRMED | closed form .822179 = integrate .822179 | - |
| 1c | "leading-order" cell 5 P_hat=.663 (bias -.337), cell 6 .238 | UNSUPPORTED as labeled | The formula is a ratio-of-means heuristic (E sqrt|X2| plugged inside 1/(a+b)), not a leading-order expansion. Literal leading term -c*sqrt(s2)/sqrt(d1) is -.508 (cell 5) and +.313 (cell 6). Quadrature bias at cell-5 s: -.325 (P_hat .675); heuristic .663 is close by luck (b/a=.51, not small). bias/leading term is .64 at s=.0428, reaching .964 only at s/256. Cell 6 quadrature used s1/s2=.522 but archived ratio is .432 (cell 6 s1=.01726, s2=.0399); heuristic .238 vs quad .235, harmless | low |
| 1d | A1 slope "4 smallest s" = .459 | MISLABELED | script uses `last4 <- 5:8` (not 6:9; s2g has 9 entries). 6:9 gives .470. Still passes; local slopes .349,.377,...,.479 rise toward .5. A2 uses index 9, ratio .964 reproduced | low |
| 1e | zero-contrast bias is O(s^(1/2)) | CONFIRMED | slopes above; sqrt non-differentiable at 0 | - |
| 2a | snap formula (independence product) and exact bivariate-normal version | CONFIRMED (as a fixed-V Gaussian statement) | cell 2 mu=(5.434,2.729): pw2*(1-pw1)=.2209 vs archived .2295. Formula xor = pw1+pw2-2*p11 is correct; independence gap <.004 per T2b. Ignores studentization (v_hat random, cor(d,v_hat)=.72-.77), which is why archived .994/.997 exceed formula .983 | low |
| 2b | P_tilde = k4 + (c-k4)S; bias = bias(k4)+E[(c-k4)S] | CONFIRMED but vacuous | Algebraic identity, no content; c depends on which axis snaps and S is dependent on k4, so E[(c-k4)S] does not factor. "Decomposes the bias" is bookkeeping, not a derivation | low |
| 2c | cell 2 "+.096 bias because snap sends estimate to 1 when truth .51" | LOOSE | .22*(1-.5097)=.108, archived .0964 (rest from snap-to-0 and non-snap bias) | low |
| 3.1 | cell 5: gate => E5 >= .9901 | CONFIRMED | (1-E)<.1 sqrt(E(1-E)) <=> E>1/1.01=.990099. Strict ">" not ">=" (immaterial). SD<=sqrt(E(1-E)) valid since Var=E T^2-E^2<=E-E^2 for T in [0,1]. Edge: SD=0 gives 0/0, gate returns NA/fail | - |
| 3.2 | cell 7: E7 <= .811 (loose .8221) | CONFIRMED | root of x-.7721=.1 sqrt(x(1-x)) is .81124; loose .7721+.05=.8221. Lower-side (E7<.7721) irrelevant for the needed TV | - |
| 3.3 | |E5-E7|<=TV; needed TV>=.179 | CONFIRMED | .990099-.81124=.17886 (strict). Valid for randomized T in [0,1] on the full data law | - |
| 4a | LR formula TV=E_P[(1-Q/P)_+] as coded | CONFIRMED | correct; LLR sums Bernoulli terms with Q (shift) in numerator, sampled under P=cell 5; E_P[Q/P]=.9992 sanity | - |
| 4b | cells 5 and 7 differ only by .05*A in the Y logit | CONFIRMED | validate script: alpha .7 both, interaction 0 both, direct 0 vs .05, n=500 both; X=(C,A,M) law identical | - |
| 4c | TV=.1404 | CONFIRMED | independent run (seed 20261021, R=40000): E_P (1-Q/P)+ = .1407+-.0008; E_Q (1-P/Q)+ = .1412+-.0008; Gaussian-LLR approx .1413. Well below .179 (>40 SE) | - |
| 5a | s^2 <= R/(R-1) m(1-m) for values in [0,1] | CONFIRMED | sum(x-m)^2=sum x^2-Rm^2<=Rm(1-m) | - |
| 5b | "worst-case SD .0114" and Bernstein <=.0042 | WRONG derivation, conclusion true | SD .0114 plugs E5=.9901 into Var5<=E5(1-E5); that is not worst case over E5 (E5=.75,E7=.61 gives difference-SD .0146 and needed t .039 = 2.7 SD). The difference-only event is the wrong object: pass needs both m5>=.99009 and m7<=.81124 individually. Using the product (independent cells), worst case over E5 in [.80,.99] gives normal <=6.4e-7, Bernstein <=8.3e-5, so P(pass)<=.0042 is true but loose and the stated 3.4 SD argument is invalid as worst-case | med |
| 5c | independence of cells | CONFIRMED | per-cell seeds 20461020+cell*10000+r distinct | - |
| 6a | psi homogeneous deg 2 in (d,s); P_hat scale invariant; axis symmetry | CONFIRMED (jointly scaling d and s; swap needs v swapped too) | a, tau*v both scale as lambda^2 | low |
| 6b | kappa=1 makes d^2-v unbiased for delta^2 | UNSUPPORTED | holds only for fixed/independent v; own notes: cor(d_hat1,v_hat11)=.72-.77, so v_hat is not independent | low |
| 6c | zero-contrast coefficient E[m]/sqrt(s) | CONFIRMED | change of variable d=s z; kappa=0,tau->0 recovers .822 | - |
| 6d | "at most 2x" coefficient, "1.7x" cells 5/6, "tenfold" regular | MOSTLY CONFIRMED | .940->.420 = 2.24x (vs raw .822: 1.96x). Cell 5 -.365->-.213 = 1.71x, but cell 6 .268->.186 = 1.44x. "Tenfold" holds vs the kappa=0 row (11.3x cell 2, 8.1x cell 4), but vs k4 (stated baseline) it is 18x and 28x. Baseline inconsistent | low |
| 6e | b/SD ranges (2-11 cells 5,6; .22 cell 7; regular fail kappa>=1) | CONFIRMED vs output file | cells 5/6 range 1.98-10.73; cell 7 best .217 (sign crossing, bias +.013); kappa=0,tau=.5 already fails cell 2 (.103) | - |
| 7-Q1 | two-point minimax floor from DGP laws | VALID concept (Le Cam) | needs two data laws; floor over a lone cell is 0 (constant=truth) and is the right reason a single cell cannot be gated by risk | - |
| 7-Q2 | max risk over ladder of weak signals | PARTLY VALID | blocks Hodges-type snap, but ladder is finite so risk between points unbounded (Codex G finding right). Ladder "in SE units" of d_hat2 understates distinguishability: TV .1404 equals Delta=.353 sd units of a Gaussian shift (2*Phi(Delta/2)-1), vs .2285 SE of d2; ratio ~1.55 (Cauchy-Schwarz, path not least favorable). Floors must come from data-law TV along the actual DGP path, not SE units | med |
| 7-Q3 | tolerance chosen on scientific grounds then checked against floor | VALID logic; unsupported until floor computed | no math content | - |
| 7-Q5 | RMSE scales as n^(-1/4) | OVERSTATED (reversal justified) | D1 is a bias rate at an exact zero; at regular cells RMSE is n^(-1/2). n^(-1/4) is correct only as a minimax/local-ladder rate (two-point floor gap ~ sqrt(k s2/d1), d1 fixed) | med |
| 7-Q6 | withheld SE near flag | no math claim; selective-testing concern stands | - | low |
| 7-D3 | "D3 is argument, not theorem" | UNDERCLAIMED | two-point (Le Cam) gives a theorem: P(d1,0)=1 vs P(d1,s)=1-O(s^(1/2)) with TV bounded away from 1 => any estimator (any threshold, even lambda_n) has sup risk >= c n^(-1/4). Needs only the Gaussian limit experiment | med |

## Missed by the authors
- Bounded-range caveat: gate bias/SD is passable by any estimator outside [0,1] via added noise; the infeasibility is a property of the bounded-range, scale-free gate. State that explicitly (one line).
- The T2b A1 index bug (5:8 vs 6:9) and cell 6 ratio .522 vs .432 in the rate script.
- Cell 7 smooth-candidate "bias +.013" is a sign crossing from over-shrinkage, not accuracy (b/SD .217 still fails).
- "kappa=1 unbiased" and the Gaussian-fixed-V snap formula both ignore the v_hat dependence the authors themselves document.
- Fix lines: redo the empirical bound with the product of the two marginal events; compute floors from data-law TV along the DGP path; promote D3 to a two-point lemma.
