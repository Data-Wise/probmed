# MBCO: which optimization method, and what medfit does (2026-10-14)

## What medfit uses

medfit's native SEM engine (`R/sem-optim.R:304`, v0.10.0) uses **`NLOPT_LD_SLSQP`
with an analytic gradient**, with `xtol_rel = 1e-10`, `ftol_rel = 1e-12` and
`maxeval = 10000`. Around it, medfit adds:
- preconditioning;
- perturbed multistarts, keeping the best accepted start when the model has
  a nonlinear constraint;
- a KKT acceptance gate, which checks for negative variances and the gradient
  sign at an active bound;
- a polish step.

SLSQP handles the constraints as general nonlinear constraints. Without
variance bounds, that lets it reach negative variances (Data-Wise/medfit#213).

## Ways to compute a likelihood-ratio (MBCO) interval

| Method | How endpoints are found | Strengths | Weaknesses for P_med |
|---|---|---|---|
| **Profile + inversion** (current probmed) | Profile LR on a grid, then `uniroot` | Detects disconnected acceptance sets; handles the composite null {a=0}∪{b=0} exactly; verified statuses | Many profile evaluations, so slower |
| **Neale & Miller (1997)** (OpenMx `mxCI`) | Optimize the parameter subject to 2ΔLL = χ²; SLSQP/NPSOL | One constrained solve per endpoint; parametrization-invariant | Assumes a connected interval, and can return one piece of a disconnected set; weak near boundaries without the Wu & Neale (2012) / Pritikin et al. (2017) correction; the kink at the composite null |
| **Venzon & Moolgavkar (1988)** | Modified Newton on the endpoint equations | Very fast when the likelihood is locally quadratic | Fragile when non-convex or poorly identified; not directly for functions of parameters (P_med is one); a 2021 trust-region fix (RVM) exists |
| **General nonlinear-constraint solve** (medfit `fit_sem` + `RMediation::mbco`) | SLSQP fit for each constrained model, compared by LR | Reuses general SEM machinery | Needs variance bounds (medfit#213); still needs an outer inversion |

## What makes probmed's problem easy

After the constraint is eliminated exactly, the inner problem is a
**smooth 2-variable unconstrained maximization** over (log Vm, log r). `a`
comes from the constraint, `b` is the closed-form profile root, and the linear
nuisance terms are profiled by QR. For a*b the inner problem has one variable.
Any reasonable optimizer reaches the same optimum: `optim` and nloptr agree to
1e-11 on 84 fits.

The analytic gradient is also cheap. Because `b` is profiled
(∂ℓ/∂b = 0), the envelope theorem gives ∂ℓ/∂(log r) from the partial
derivative only, and ∂ℓ/∂(log Vm) follows directly through
a(Vm, r) ∝ √(Vm + r).

## Recommendation

1. **Keep profile + verified inversion as the outer method.** It is the only
   option here that handles the composite null and detects disconnected
   acceptance sets.
2. **Inner solver: match medfit.** Use `NLOPT_LD_SLSQP` (or `NLOPT_LD_LBFGS`,
   since the reduced problem has no constraints) with the analytic gradient of
   the 2-variable profile, and keep BOBYQA as the fallback. This should more
   than recover the 1.7× slowdown in #48 and puts both packages on one
   gradient-based nloptr family.
3. **Add an independent oracle to the tests.** Neale–Miller endpoint search
   via `RMediation::mbco` (or OpenMx) on interior cases, as an
   implementation-independent check (Suggests only).
4. **Parallel-mediator MBCO** uses a penalty wall instead of exact
   elimination, and is the remaining weak spot. Apply the same elimination
   there, or bound the variances and solve it as a constrained problem with
   SLSQP, as medfit does.

## Sources
- Neale, M. C., & Miller, M. B. (1997). The use of likelihood-based confidence intervals in genetic models. *Behavior Genetics, 27*, 113–120. OpenMx `mxCI`: https://search.r-project.org/CRAN/refmans/OpenMx/html/mxCI.html
- Pritikin, J. N., Rappaport, L. M., & Neale, M. C. (2017). Likelihood-based confidence intervals for a parameter with an upper or lower bound. *Structural Equation Modeling, 24*(3), 395–401. https://pmc.ncbi.nlm.nih.gov/articles/PMC5863926
- Venzon, D. J., & Moolgavkar, S. H. (1988). A method for computing profile-likelihood-based confidence intervals. *JRSS C, 37*(1), 87–94. https://ideas.repec.org/a/bla/jorssc/v37y1988i1p87-94.html ; RVM extension: https://link.springer.com/article/10.1007%2Fs11222-021-10012-y
- Tofighi, D. (2020). Bootstrap model-based constrained optimization tests of indirect effects. *Frontiers in Psychology*. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6984355/
