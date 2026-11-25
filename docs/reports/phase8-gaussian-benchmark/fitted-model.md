# Local calibration of the fitted Gaussian profile LR

## Result and scope

Consider independent homoskedastic Gaussian errors in the single-mediator
models

$$
M=W_n\alpha+\varepsilon_m,\qquad
Y=W_n\gamma+bM+\varepsilon_y.
$$

Condition on a full-rank design \(W_n\) of fixed column count \(q\) with
treatment \(x\) and the same controls in both equations. Write \(a\) for
the treatment coefficient in \(\alpha\). Fix positive true residual
variances \(V_{m0},V_{y0}\) and a treatment contrast \(\Delta\ne0\).
Let \(A_n\) be the squared norm of treatment residualized on the other
columns of \(W_n\), and assume \(A_n/n\to d>0\).

The modeled probability is

$$
p=\Phi\!\left(\frac{\Delta ab}{\sqrt{2(b^2V_m+V_y)}}\right).
$$

Write \(L_n(p)\) for twice the difference between the unrestricted and
globally constrained maximized log likelihoods. For any fixed compact
path region \(K\subset\mathbb R^2\),

$$
\limsup_{n\to\infty}\sup_{(\sqrt n\,a,\sqrt n\,b)\in K}
\Pr\!\left[L_n(p_{\rm true})>\chi^2_{1,.95}\mid W_n\right]\le .05.
$$

Consequently the exact inversion has asymptotic coverage of at least
95% over such regions, with coverage tending to 99.75% at \(a=b=0\).
The model permits common covariate adjustment and fits both residual
variances. This argument concerns exact likelihood maximization and
inversion. The package uses numerical profiling and an endpoint scan.
[Gaussian benchmark](proof.md) · [Profile code](../../../R/compute-mbco-verified.R)

## Profile likelihood and Gaussian coordinates

Let \(B_n=\|(I-P_{W_n})M\|^2\), and write \(\widehat V_m,\widehat V_y\)
for the unrestricted residual sums of squares divided by \(n\). After
profiling the linear nuisance coefficients, twice the log-likelihood
loss at \((a,b,V_m,V_y)\) is exactly

$$
D_n=\sum_{j=m,y} n\left\{\log(V_j/\widehat V_j)
       +\widehat V_j/V_j-1\right\}
   +\frac{A_n(a-\widehat a)^2}{V_m}
   +\frac{B_n(b-\widehat b)^2}{V_y}.
$$

This follows by completing the square in each residual sum of squares.
In the mediator regression, \(B_n=n\widehat V_m\) and
\(B_n/V_{m0}\sim\chi^2_{n-q}\). In the outcome regression, the standardized
\(b\) error
\(Z_b=\sqrt{B_n}(\widehat b-b)/\sqrt{V_{y0}}\) is standard normal
conditional on \(M,W_n\). Its conditional law is constant, so \(Z_b\)
is independent of the mediator sample conditional on \(W_n\). The
standardized \(a\) error is standard normal by fixed-design Gaussian
regression. Both ML variance estimates are root-\(n\) consistent, by
the residual chi-square laws with \(n-q\) and \(n-q-1\) degrees of freedom.
[Profile setup](../../../R/compute-mbco.R)

Set \(s_a=\sqrt{V_{m0}/d}\), \(s_b=\sqrt{V_{y0}/V_{m0}}\), and define

$$
u=\sqrt n\,a/s_a,\quad v=\sqrt n\,b/s_b,\quad
w=v\sqrt{\frac{V_{y0}}{b^2V_m+V_y}}.
$$

In \((u,w)\) the constraint at a target \(p_n\) is exactly

$$
uw=\theta_n=\frac{n\sqrt{2V_{y0}}\,\Phi^{-1}(p_n)}
                   {\Delta s_a s_b}.
$$

For bounded \(w\) and variances near their positive true values, the
inverse is \(b=t\sqrt{V_y}/\sqrt{V_{y0}-t^2 V_m}\), where
\(t=s_b w/\sqrt n\). This inverse also works for \(w=0\). Thus both
branches at every nonzero target and both axes at the zero target are
represented by the exact constraint.

## Transfer of the size bound

Let \(\lambda_j=\sqrt n(V_j/\widehat V_j-1)\) and
\((\widehat u,\widehat v)=(\sqrt n\,\widehat a/s_a,
\sqrt n\,\widehat b/s_b)\). On bounded sets in \((u,w,\lambda)\),
\(w-v=o_p(1)\) uniformly because the fitted variances are consistent
and \(b^2=O(n^{-1})\). Expansion of the variance losses and
\(A_n/n\to d,\ B_n/n\to V_{m0}\) give

$$
D_n=(u-\widehat u)^2+(w-\widehat v)^2
       +\tfrac12(\lambda_m^2+\lambda_y^2)+o_p(1)
$$

uniformly on these sets. The expansion uses
\(\log r+r^{-1}-1=(r-1)^2/2+O(|r-1|^3)\) near \(r=1\).

At \(p_n=p_{\rm true}\), the true parameter satisfies the constraint
and has \(D_n=O_p(1)\) uniformly on the stated local regions. Each
term of \(D_n\) is nonnegative. The variance loss is zero only at
\(r=1\) and diverges at both ends of \(r>0\). Thus points with
bounded loss have \(V_j/\widehat V_j=1+O_p(n^{-1/2})\).
The path losses then bound both standardized path coordinates in
probability. It follows that constrained near minimizers may be
localized to bounded \((u,w,\lambda)\) sets with arbitrarily high
probability. The quadratic expansion and the exact constraint give
a lower bound by the squared distance to \(H_{\theta_n}=\{uw=\theta_n\}\)
up to \(o_p(1)\).

For the upper bound, choose the closest point of \(H_{\theta_n}\) to
\((\widehat u,\widehat v)\). It is bounded in probability since
the true standardized point is feasible and bounded. At the closest
point set \(V_j=\widehat V_j\), and use the inverse map to produce
a feasible \((a,b)\). With \(\lambda_j=0\) the loss achieves the
squared distance up to \(o_p(1)\). Thus

$$
L_n(p_{\rm true})=
\operatorname{dist}^2((\widehat u,\widehat v),H_{\theta_n})+o_p(1).
$$

The estimates equal \((u_{\rm true},w_{\rm true})+Z+o_p(1)\)
with \(Z\sim N_2(0,I)\) by the regression representation above.
These probability bounds are uniform over bounded true local
paths. For example, \(B_n/n\), the two standardized errors, and
the two ML variance ratios have laws independent of \(a,b\) under
the common design. Also \(v_{\rm true}-w_{\rm true}=O(n^{-1})\)
uniformly in that region. Distance to a fixed closed set is
1-Lipschitz in the observation. So the same LR approximation holds
with \((\widehat u,\widehat v)\) replaced by the Gaussian observation.

For the rejection probability, along any subsequence the bounded
true Gaussian means have a convergent subsequence, and so do their
product targets. The distances are continuous in the observation
and target. For the target coordinate, a limit of bounded closest
points is feasible at the limit target. Conversely, a point
\((s,t)\) of the limit curve with \(s\ne0\) is approximated by
\((s,\theta_n/s)\), and there is a symmetric construction when
\(t\ne0\). At the origin use
\((\sqrt{|\theta_n|},\operatorname{sgn}(\theta_n)\sqrt{|\theta_n|})\).
This supplies both distance bounds, including when the target
approaches zero.

There is no mass at the positive critical distance. For a nonzero
target, points at that distance lie on the fixed-radius normal
offsets of the smooth hyperbola branches. These curves have
Lebesgue measure zero. At target zero the level set lies on lines
parallel to the axes. The Gaussian benchmark therefore bounds
the limiting rejection by 0.05 along every such subsequence,
which proves the uniform bound. At the double null the limiting
statistic is the smaller of two independent squared standard
normals. This gives limiting rejection \(0.05^2\).

## Remaining work

- Establish control for sequences leaving every bounded local path
  region, including transition cases where one path strengthens.
- Check convergence of the constrained numerical maxima and inversion
  endpoints on both branches at the local \(1/n\) target scale.
- A smaller interval still requires a data-dependent test family and
  the frozen lattice checks. Phase 8 remains NO-GO.
