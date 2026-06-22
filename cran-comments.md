## Release summary

probmed 0.3.0 is distributed as a **GitHub (non-CRAN) release**. It adds four
estimator families for modern proportion-mediated effect sizes (gauge-calibrated
`ward_residual()`, incremental `incr_pmed()`, Sobol `sobol_pmed()`, and the shared
`pmed_sensitivity()` helper), all built on a cross-fitted corner-EIF core.

## R CMD check results

0 errors | 0 warnings | 1 note

The single NOTE is expected and intentional for this GitHub release:

```
checking CRAN incoming feasibility ... NOTE
  New submission
  Unknown, possibly misspelled, fields in DESCRIPTION:
    'Remotes'
```

* **`Remotes`** — `Remotes: data-wise/medfit@v0.3.0` pins the foundation package
  `medfit`, which is not yet on CRAN. The field is load-bearing: it lets
  `remotes::install_github("data-wise/probmed")` resolve `medfit` at the correct
  version. It is **not** a CRAN-recognized field, hence the NOTE. **This field
  will be removed at the time of CRAN submission**, once `medfit` is itself on
  CRAN and can be declared as an ordinary `Imports` dependency.
* **New submission** — informational; the package is not currently on CRAN.

No further action is required for the GitHub release. The CRAN submission is
gated on `medfit` reaching CRAN (tracked in the mediationverse ecosystem plan).
