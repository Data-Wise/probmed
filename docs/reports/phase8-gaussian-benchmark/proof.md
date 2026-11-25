# Sharp 95% profile-likelihood calibration for a product of Gaussian means

## Result

Let $X=(X_1,X_2)\sim N_2(\mu,I_2)$, with $\mu=(u,v)$, and write
$H_\theta=\{(a,b):ab=\theta\}$ and
$T_\theta(X)=\operatorname{dist}(X,H_\theta)^2$. This is twice the
log-likelihood difference between the unrestricted and constrained fits.
Let $k=\Phi^{-1}(.975)$. For each real $\theta$,
$$
\sup_{\mu\in H_\theta}\Pr_\mu(T_\theta(X)>k^2)=.05.
$$
Consequently $k^2=3.84145882\ldots$ is the smallest constant critical value
that controls size for each product target. The proof uses geometry outside a
disk of mean values and finite exact arithmetic inside it. It applies to the
known identity-covariance Gaussian experiment.

For $\theta=0$, $T_0(X)=\min(X_1^2,X_2^2)$, and at least one of the
coordinates has mean zero. This gives the upper bound directly. For
$\theta\ne0$, coordinate sign changes and swapping preserve all distances and
reduce the proof to $u\ge v>0$ and $\theta>0$.

## Means with norm at least $k/\sqrt{2}$

Fix such a mean and put $m=\|\mu\|$. Conditional on the length
$r=\|X-\mu\|$, the error direction is uniform. When $r\le k$, the point
$\mu\in H_\theta$ alone gives acceptance. When $r>k$, put
$$
L=\sqrt{r^2-k^2},\qquad \gamma=\arcsin(k/r).
$$
On each side of $\mu$ on its hyperbola branch there is a point at distance
$L$ from $\mu$. This follows from continuity and the divergence of norm
at both ends of the branch, parameterized by $(a,\theta/a)$ for $a>0$.
A closed $k$-disk around either point covers an arc of half-angle $\gamma$
on the radius-$r$ observation circle. The cosine law gives the cut
$\cos(\delta)\ge L/r$ for angular deviation $\delta$ from that chord.

The chords have angular separation at least $2\gamma$, as follows.

- If $r\ge\sqrt{2}k$, then $2\gamma\le\pi/2$. The two chords lie in
  opposite coordinate quadrants by monotonicity of $\theta/a$, so their
  angular separation is at least $\pi/2$.
- If $k<r<\sqrt{2}k$, use orthonormal coordinates
  $S=(a+b)/\sqrt2,\ W=(a-b)/\sqrt2$, in which the branch is
  $S=\sqrt{2\theta+W^2}$. Write the rotated mean as $(s,w)$. The unit
  normal is $n=(s,-w)/m$ and the unit tangent is $\tau=(w,s)/m$.
  For a chord $d=(d_S,d_W)$ of length $L$, subtraction of the two
  hyperbola equations gives
  $2m\,n\cdot d=d_W^2-d_S^2$, so
  $|n\cdot d|/L\le L/(2m)$.
  The branch graph has slope of absolute value less than one, so
  $|d_S|<|d_W|$. Since $s>|w|$, the sign of $\tau\cdot d$ is the sign
  of $d_W$. Hence the two chords are within
  $\arcsin(L/(2m))$ of opposite tangent directions.
  Here $r\le2m$, so their separation is at least
  $\pi-2\arcsin(L/(2m))\ge\pi-2\arcsin(L/r)=2\gamma$.

Thus the two arcs have disjoint interiors and cover $4\gamma$ radians.
That equals the angular length within distance $k$ of any straight line
through $\mu$ on this same observation circle. Averaging over $r$ shows
that the hyperbola tube has probability at least that of the line tube,
which is $2\Phi(k)-1=.95$.

## Means with norm less than $k/\sqrt{2}$

The certificate covers all $u\ge v\ge0$ with $u^2+v^2\le 2$, a region
containing the remaining means since $k<2$. It uses a smaller distance
threshold, $1.95<k$. The two bounds on $k$ follow from alternating Taylor
bounds for the normal integral. Specifically set
$$
B_n(x)=\sum_{j=0}^{n}\frac{(-1)^j x^{2j+1}}{2^j j!(2j+1)}.
$$
Taylor's remainder for `exp(-t^2/2)` has sign $(-1)^{n+1}$ for
$t\ge0$, giving upper and lower integral bounds at even and odd $n$.
Using $0.398<(2\pi)^{-1/2}<0.399$, which follows from
$223/71<\pi<22/7$, exact arithmetic gives
$$
\tfrac12+0.399 B_{18}(1.95)<.975<\tfrac12+0.398 B_{19}(2).
$$

Partition the positive mean quadrant into closed boxes of side $0.1$,
retaining the 92 boxes with lower corner $(u_0,v_0)$ satisfying
$u_0\ge v_0\ge0$ and $u_0^2+v_0^2\le2$. This covers the required
parameter region. Partition each error coordinate into 120 bins of length
$0.1$ on $[-6,6]$. Coordinates in the certificate are scaled by 100.

At an observation $X=\mu+z$, consider these 17 displacements $d$:
the origin; the four axis points at distance 1.95; the four independent
signs of $(1.37,1.37)$; and the independent signs of $(1.80,0.75)$
and its swap. All have norm at most 1.95. With $a=z_1+d_1$ and
$b=z_2+d_2$, the residual product is
$$
g_d(u,v,z)=(X_1+d_1)(X_2+d_2)-uv=ab+ub+va .
$$
It is affine separately in each of the four arguments. On each product of
a mean box and an error box, look for one $d$ whose maximum residual is
nonpositive and one whose minimum is nonnegative. Both extrema occur at
vertices. If both witnesses exist, the segment between the displaced
observations stays inside the ball of radius 1.95 around $X$ and contains
a point of $H_{uv}$, by continuity of the product. That entire parameter
and error box therefore accepts at critical radius $k$.

For each error bin $J$, let $e$ be its nearest endpoint to zero and
$D(x)=\sum_{j=0}^{20}x^j/j!$. Its Gaussian probability is bounded above
by $0.1\cdot0.399 / D(e^2/2)$. This follows from monotonicity of the
density away from zero and $D(x)\le\exp(x)$ for $x\ge0$. Round this
probability up to a multiple of $10^{-8}$ and call the resulting integer
numerator $w(J)$. The two-dimensional box mass is bounded above by the
product of the two weights divided by $10^{16}$, by independence. Mills'
inequality and the union bound control the probability outside the error
square by $4\cdot0.399/[6D(18)]$. Its upper integer numerator on the
same scale is $55{,}440{,}850$.

For each parameter box, add this tail numerator and the products
$w(J_1)w(J_2)$ over all error boxes without both witnesses. The largest
sum is $264{,}990{,}180{,}241{,}388$ in the parameter box with lower corner
$(1.3,0.5)$. Dividing by $10^{16}$ gives an upper bound of
$0.0264990180241388<.05$. This certifies every parameter in the retained
boxes. The sum uses integer and rational arithmetic, with no root finding,
floating-point decisions, or simulation. The companion Python certificate
computes extrema by choosing the mean endpoints according to signs. A
separate C++ implementation checked all 16 vertices directly and reproduced
the maximum numerator using the same exact bin weights.

## Sharpness and inversion

For any fixed $\theta$, take $\mu_m=(m,\theta/m)$ with $m\to\infty$
and couple observations with one fixed standard normal error $z$. For large
$m$, the point $(m+z_1,\theta/(m+z_1))$ gives
$\limsup T_\theta(\mu_m+z)\le z_2^2$. A nearest point in $H_\theta$
exists because that set is closed. Its first coordinate is $m+O(1)$ because
the distance just bounded is $O(1)$, so its second coordinate tends to zero.
The squared vertical difference thus gives
$\liminf T_\theta(\mu_m+z)\ge z_2^2$. This proves convergence to
$\chi_1^2$, including when $\theta=0$. Its continuity at the critical
value proves sharpness of the uniform size.

The inverted set is the image of the closed disk of radius $k$ about $X$
under $(a,b)\mapsto ab$. This image is a compact interval by continuity
and connectedness, with infimum coverage 95% in this experiment. Extending the
claim to a fitted mediation model requires a separate link from its profiled
likelihood to this Gaussian experiment, uniform over the relevant local
parameter sequences.

## Computational checks

Run `python3 certify_weak.py` to reproduce the exact certificate. The
separate `grid.R` calculation uses deterministic numerical integration of
vertical tube slices and is an exploratory check only. Its 3,200 pointwise
evaluations found no exceedance of .05.
