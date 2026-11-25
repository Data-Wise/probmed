# Gaussian product test family with a radial rejection rule

**2026-11-25.** This pilot extends the fixed-target rule to every real
product in the independent Gaussian shift experiment. It gives at least
95% coverage. The origin and strong-path results below explain why this
rule stays a benchmark for Phase 8 research.

## Rule and uniform size

Let $X\sim N_2(\mu,I_2)$, with product target $\theta=\mu_1\mu_2$.
For each hypothesized $\theta$, define $H_\theta=\{(s,t):st=\theta\}$,
$T_\theta(X)=\operatorname{dist}^2(X,H_\theta)$, and $S=\|X\|^2$.
Choose $0<\eta<0.05$. For $q(c)=\Pr(\chi^2_1>c)$, define

$$
G(c)=\frac{q(c)^2+2q(c)+e^{-c/2}}4,qquad G(c)=0.05-\eta.
$$

Use this single $c$ for all targets. Let $d(\theta)$ be the $\eta$
quantile of $\chi^2_2$ with noncentrality $2|\theta|$. Reject $\theta$ if
$T_\theta(X)>c$ or $S<d(\theta)$.

**Distance bound.** For $\theta>0$ and $\mu\in H_\theta$, both
directions along the branch through $\mu$ continue to arbitrarily large
distances. Monotonicity of $t=\theta/s$ puts the two chords from $\mu$
in opposite coordinate quadrants, with angle $\beta\in[\pi/2,\pi]$.
Reflection gives the same geometry for $\theta<0$. At $\theta=0$ choose
the two directions on an axis containing $\mu$.

Condition on $r=\|X-\mu\|$. If $r\le\sqrt c$, all directions are
accepted by the distance rule. Otherwise choose one point on each side
at chord length $\sqrt{r^2-c}$. Their radius-$\sqrt c$ disks each cover
an arc of half-angle $\gamma=\arcsin(\sqrt c/r)$. The two arcs cover
length $4\gamma-\max(0,2\gamma-\beta)$, which is least when
$\beta=\pi/2$. Those two arcs at right angles give exactly the
directions within distance $\sqrt c$ of two perpendicular rays from
$\mu$. For the rays, the squared distance has tail $q^2$ in their
shared quadrant, $q$ in each adjacent quadrant, and $e^{-c/2}$ in the
opposite quadrant. Signs are independent of squared magnitudes, so their
equal mixture gives $\Pr_\mu(T_\theta>c)\le G(c)$. Comparing the same
arcs with opposite rays shows $G(c)\ge q(c)$, since opposite rays form
a line. In particular $c>\chi^2_{1,0.95}$.

**Radial bound.** On $H_\theta$, $\|\mu\|^2\ge 2|\theta|$.
The Gaussian mass of an origin-centered disk decreases with $\|\mu\|$.
Rotate to mean $(m,0)$. Given the second coordinate, the disk slice is
empty or $[-a,a]$, whose normal probability has derivative
$\phi(a+m)-\phi(a-m)\le0$ for $m\ge0$. For a disk of positive radius,
the integrated probability decreases strictly for $m>0$. Therefore
$\Pr_\mu(S<d(\theta))\le\eta$. The union bound gives uniform rejection
probability at most $G(c)+\eta=0.05$ on each $H_\theta$.

## Inversion and empty sets

The distance condition accepts precisely the product image of the
closed disk of radius $\sqrt c$ about $X$. A continuous function on this
compact connected disk has a compact interval as its image. The radial
condition accepts a symmetric closed interval when $S\ge d(0)$ and
the empty set otherwise. Indeed $d(\theta)$ is continuous, strictly
increasing in $|\theta|$, and diverges as $|\theta|\to\infty$. These
quantile properties follow from the disk probability and the continuous,
strictly increasing CDF in the disk's squared radius. Intersecting the
two sets gives a closed interval whenever the result is nonempty.

In fact, the inversion is empty **exactly when $S<d(0)$**. For
$S\ge d(0)$ and $\min(|X_1|,|X_2|)\le\sqrt c$, the two conditions
accept zero. Otherwise $S>2c>4$. Take $\theta=X_1X_2$, so
$T_\theta(X)=0$ and $2|\theta|\le S$. For $m=\sqrt S>2$ and
$Z\sim N_2(0,I_2)$, the event $-1\le Z_1\le-1/2,\ |Z_2|\le1$
implies $\|Z+(m,0)\|^2\le S$ and has probability greater than 0.10.
By disk monotonicity the CDF at $S$ with noncentrality $2|\theta|$ is
at least this large, so $d(\theta)\le S$. Thus $\theta=X_1X_2$
is accepted in this case. At the origin, the exact empty-set probability
is $\eta$.

## Origin coverage and strong-path cost

At $\mu=0$, the distance rejection event at $\theta=0$ requires both
$|X_i|>\sqrt c$. It is disjoint from the disk rejection event because
$d(0)<0.103<c$. Origin coverage is exactly $1-\eta-q(c)^2$. It is
99.75% for the ordinary $\chi^2_1$ inversion.

For a strong-path sequence $\mu=(m,v_m)$ with $m\to\infty$ and bounded
$v_m$, couple all observations as $(m+Z_1,v_m+Z_2)$. The true target is
$mv_m$. Evaluating the distance at $(m+Z_1,mv_m/(m+Z_1))$ gives an
upper bound approaching $Z_2^2$. A closest point has first coordinate
$m+O(1)$ since the distance to the true mean is $\|Z\|$. Its second
coordinate is $v_m+O(1/m)$, giving the same lower limit. Signs and
coordinate swaps cover either strong path. The radial test is eventually
inactive at the true target. By Markov's inequality
$d(\theta)\le(2+2|\theta|)/(1-\eta)$, while $S$ grows as $m^2$.
Coverage thus tends to $1-q(c)$ along these sequences.

The same radial bound makes the radial condition eventually inactive
throughout the distance interval, where $|\theta|=O(m)$. For a disk of
radius $k$, the product at displacement $h$ from $X$ equals
$X_1X_2+X_1h_2+X_2h_1+h_1h_2$. The last two terms are bounded by
$k|X_2|+k^2/2$, and the points $(h_1,h_2)=(0,\pm k)$ attain both linear
extremes $\pm kX_1$. Hence the interval width is $2k|X_1|+O(1)$.
Its ratio to the $\chi^2_1$ width tends to $\sqrt{c/\chi^2_{1,0.95}}$.

| $\eta$ | $c$ | Origin coverage | Empty at origin | Strong-path limit of coverage | Strong-path width ratio |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.010000 | 4.659945 | 98.905% | 1.000% | 96.913% | 1.1014 |
| 0.029782 | 5.945960 | 97.000% | 2.978% | 98.525% | 1.2441 |
| 0.040000 | 7.291795 | 95.995% | 4.000% | 99.307% | 1.3777 |

The middle row is the minimum $\eta$ in this family that brings ideal
origin coverage down to 97%. To see this, write the coverage as
$0.95+G(c)-q(c)^2$. It decreases in $c$ because its derivative is
$(1/2-3q(c)/2)q'(c)-e^{-c/2}/8<0$ for $q(c)<0.05$.
Increasing $\eta$ increases $c$, the empty-set probability, and the
limiting width ratio. Meeting the 97% ideal-origin target therefore
entails at least 2.978% empty sets there and a 24.41% limiting width
increase for this construction. The exact ideal-origin probability is
separate from the frozen T5 per-cell screening bound.

## Deterministic checks and disposition

Run `Rscript --vanilla check_test_family.R` from this directory.
The script integrates the vertical slices of the distance tube, with
the origin disk removed for the pilot. Branch extrema are obtained from
the slice boundaries, the limit at $s=0$, and the stationary points of
$\theta/s\pm\sqrt{c-(s-x)^2}$. The twelve fixtures cover sign and
coordinate symmetries, axes, balanced and unequal paths, and one
off-target comparison. Two integration partitions agreed within
$9.0\times10^{-10}$. The three axis fixtures also agreed with direct
normal and noncentral chi-square probabilities within $3\times10^{-16}$.
Quadrature is a pointwise check. The uniform guarantee follows from the
size argument above.

| True mean | Tested product | $\chi^2_1$ | $\eta=0.01$ | $\eta=0.04$ | $\chi^2_2$ |
| --- | ---: | ---: | ---: | ---: | ---: |
| (0.5, 0.5) | 0.25 | 99.483% | 98.779% | 95.984% | 99.943% |
| (1, −1) | −1 | 98.774% | 98.430% | 95.943% | 99.826% |
| (2, 0.5) | 1 | 97.297% | 98.270% | 98.451% | 99.516% |
| (12, 1/12) | 1 | 95.000% | 96.913% | 99.307% | 98.562% |
| (0, 0) | 1 | 98.960% | 96.932% | 89.669% | 99.918% |

The saved T5 data lack the observed path coordinates needed to evaluate
$S$. No lattice score or fitted-model validity claim follows from the
pilot. The next research step is to use local size slack without
enlarging the strong-path cutoff, and to check whether the resulting
inversion can be kept nonempty. Phase 8 remains NO-GO. Retain the
$\chi^2_1$ interval and `near_double_null`. I7, package publishing, and
pmed submission remain on hold. No new simulations ran.
