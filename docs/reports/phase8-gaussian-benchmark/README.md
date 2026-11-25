# Gaussian benchmark for the Phase 8 interval

Archived 2026-11-25 from the reviewed `probmed-chi1.apqI1p` scratch
directory. The proof and the two certificate sources are byte-identical
to the reviewed versions. The certificate covers a continuum of means
by bounding each parameter box.

- [Sharp 95% calibration proof](proof.md) for the independent Gaussian
  shift experiment with known unit variances.
- [Fitted-model result](fitted-model.md) for the exact Gaussian profile
  LR on bounded local path regions.

## Reproduce the certificate

Run from this directory. The Python script uses the standard library,
and the C++ check requires a C++17 compiler.

```sh
python3 -B certify_weak.py
work=$(mktemp -d)
c++ -O2 -std=c++17 verify_corners.cpp -o "$work/verify_corners"
python3 -B - <<'PY' > "$work/bins.txt"
from certify_weak import BINS, SHIFTS, TAIL_WEIGHT

print(len(BINS), len(SHIFTS), TAIL_WEIGHT)
for lower, weight in BINS:
    print(lower, weight)
for first, second in SHIFTS:
    print(first, second)
PY
"$work/verify_corners" < "$work/bins.txt"
```

Both implementations give 92 parameter boxes and largest numerator
`264990180241388` on scale `10^16`. The resulting rejection upper bound
is `0.0264990180241388` in the certified region. The C++ check enumerates
all 16 vertices directly using the exact bin weights. The geometric
argument in the proof covers the remaining parameter region.

Phase 8 remains NO-GO. T5 is closed. The existing χ²₁ interval and the
`near_double_null` diagnostic remain in use. The result is a research
benchmark. Novelty review and verification of numerical inversion are
separate tasks.
