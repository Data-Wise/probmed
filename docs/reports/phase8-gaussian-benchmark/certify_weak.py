"""Exact rational upper bound for the unit-Gaussian product test near zero."""
from fractions import Fraction
from itertools import product
from math import ceil

COORD = 100
STEP = 10
END = 600
WEIGHT = 10**8
SHIFTS = [(0, 0), (195, 0), (-195, 0), (0, 195), (0, -195)]
for a, b in [(137, 137), (180, 75), (75, 180)]:
    SHIFTS.extend((sa*a, sb*b) for sa, sb in product((-1, 1), repeat=2))

# Archimedes' bounds enclose the Gaussian density normalizer.
LOWER_DENSITY_SCALE = Fraction(398, 1000)
UPPER_DENSITY_SCALE = Fraction(399, 1000)
assert 2 * Fraction(22, 7) * LOWER_DENSITY_SCALE**2 < 1
assert 2 * Fraction(223, 71) * UPPER_DENSITY_SCALE**2 > 1

def gaussian_integral_partial_sum(x, degree):
    # Taylor remainder for exp(-t*t/2) alternates in sign for t >= 0.
    term = total = x
    for j in range(1, degree + 1):
        term *= -x*x * Fraction(2*j - 1, 2*j*(2*j + 1))
        total += term
    return total

TARGET_CDF = Fraction(975, 1000)
cdf_upper_195 = Fraction(1, 2) + UPPER_DENSITY_SCALE * (
    gaussian_integral_partial_sum(Fraction(195, 100), 18)
)
cdf_lower_200 = Fraction(1, 2) + LOWER_DENSITY_SCALE * (
    gaussian_integral_partial_sum(Fraction(2), 19)
)
assert cdf_upper_195 < TARGET_CDF < cdf_lower_200
assert all(x*x+y*y <= 195**2 for x, y in SHIFTS)

def inverse_exp_upper(x):
    """Positive Taylor terms for exp(x) give a lower bound."""
    term = total = Fraction(1)
    for j in range(1, 21):
        term *= x / j
        total += term
    return 1 / total

def density_upper(at):
    return UPPER_DENSITY_SCALE * inverse_exp_upper(at*at/2)

def bin_weight(start):
    nearest_zero = min(abs(start), abs(start + STEP))
    return ceil(WEIGHT * Fraction(STEP, COORD) *
                density_upper(Fraction(nearest_zero, COORD)))

TAIL_WEIGHT = ceil(WEIGHT**2 * 4 * density_upper(Fraction(END, COORD)) /
                   Fraction(END, COORD))
BINS = [(start, bin_weight(start)) for start in range(-END, END, STEP)]

def all_accepted(z1, z2, u, v):
    """Find fixed displacements whose products bracket uv throughout a box."""
    seen_negative = seen_positive = False
    for d1, d2 in SHIFTS:
        a = z1 + d1
        b = z2 + d2
        a1 = a + STEP
        b1 = b + STEP
        # The residual is a*b + u*b + v*a, multiaffine in (a,b,u,v).
        if not seen_negative:
            seen_negative = max(
                a*b + (u + STEP*(b > 0))*b + (v + STEP*(a > 0))*a,
                a*b1 + (u + STEP*(b1 > 0))*b1 + (v + STEP*(a > 0))*a,
                a1*b + (u + STEP*(b > 0))*b + (v + STEP*(a1 > 0))*a1,
                a1*b1 + (u + STEP*(b1 > 0))*b1 + (v + STEP*(a1 > 0))*a1,
            ) <= 0
        if not seen_positive:
            seen_positive = min(
                a*b + (u + STEP*(b < 0))*b + (v + STEP*(a < 0))*a,
                a*b1 + (u + STEP*(b1 < 0))*b1 + (v + STEP*(a < 0))*a,
                a1*b + (u + STEP*(b < 0))*b + (v + STEP*(a1 < 0))*a1,
                a1*b1 + (u + STEP*(b1 < 0))*b1 + (v + STEP*(a1 < 0))*a1,
            ) >= 0
        if seen_negative and seen_positive:
            return True
    return False

def bound(u, v):
    return TAIL_WEIGHT + sum(
        w1*w2 for z1, w1 in BINS for z2, w2 in BINS
        if not all_accepted(z1, z2, u, v)
    )

if __name__ == "__main__":
    # By coordinate reflections and swapping, u >= v >= 0 suffices.
    # Include a whole box if its lower corner lies within the radius sqrt(2).
    corners = [(u, v) for u in range(0, 150, STEP) for v in range(0, u + STEP, STEP)
               if u*u + v*v <= 2*COORD**2]
    upper, corner = max((bound(*corner), corner) for corner in corners)
    print(f"Parameter boxes checked: {len(corners)}")
    print(f"Gaussian bins per axis: {len(BINS)}")
    print(f"Upper bound numerator: {upper}, denominator: {WEIGHT**2}")
    print(f"Largest bound: {upper/WEIGHT**2:.12f}")
    print(f"Lower parameter corner: {tuple(x/COORD for x in corner)}")
    print(f"Tail allowance: {TAIL_WEIGHT/WEIGHT**2:.12g}")
    assert upper*20 < WEIGHT**2
