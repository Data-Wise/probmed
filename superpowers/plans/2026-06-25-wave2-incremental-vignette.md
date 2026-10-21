# Wave 2 — incr_pmed() Vignette Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Add a probmed vignette demonstrating the incremental mediated elasticity `P_med^delta(delta)` curve, its inference, and the new `incr_sensitivity()` sensitivity path (issue #5).

**Architecture:** A single Quarto vignette `vignettes/incremental-pmed.qmd` following the established `vignettes/gauge-residual.qmd` pattern (quarto::html engine + the Windows code-only degradation guard). Self-contained DGP; renders to HTML.

**Tech Stack:** R (>= 4.1), Quarto (`quarto::html` vignette engine), knitr opts.

## Global Constraints

- R minimum 4.1.0; native pipe `|>`; explicit `package::function()` namespacing.
- Vignette front-matter must match the existing pattern exactly: `format: html`, `%\VignetteEngine{quarto::html}`, `%\VignetteEncoding{UTF-8}`, and the `requireNamespace("probmed")` → `eval = FALSE` guard chunk (Windows R CMD check child-process safety — see memory `ci-triggers-and-release-to-main`).
- Single-file deliverable; one responsibility (teach the incremental curve).
- Branch `feature/gauge-bootstrap-se` (feature branch — new files allowed).
- No pmed-modern writes.

---

## File Structure

- Create: `vignettes/incremental-pmed.qmd` — the vignette.

---

### Task 1: Write and render the incremental-pmed vignette

**Files:**
- Create: `vignettes/incremental-pmed.qmd`

**Interfaces:**
- Consumes: `incr_pmed(data.frame(A,M,Y,C), deltas, ...)` → `IncrPmedResult` (print method shows the per-δ curve); `incr_sensitivity(fit, threshold = 0)` → per-δ data frame `delta, med, tot, Pmed, tipping, tipping_threshold`.
- Produces: a rendered HTML vignette.

**Content outline (sections the vignette must cover):**
1. **Why a curve, not a number** — the incremental tilt `q_delta(C) = delta*g/(delta*g+1-g)`; `P_med^delta(delta) = med/tot` as a function of `delta`.
2. **A worked example** — a self-contained DGP with a treatment-by-mediator interaction so the curve bends; call `incr_pmed()`, print it, read the per-δ output (flat ⇒ no A×M interaction; bends ⇒ interaction present).
3. **Inference** — note the g-cross-fit one-step estimator and pointwise CIs use the full EIF including the g-score term (Kennedy 2019 Cor. 2 term II), Neyman-orthogonal in the propensity.
4. **Sensitivity** — call `incr_sensitivity(fit, threshold = 0)`; explain `tipping = -med` as the M–Y confounding bias that would zero the incremental share at each `delta`.

- [ ] **Step 1: Write the vignette file**

Create `vignettes/incremental-pmed.qmd` with this exact front-matter and guard, then the four sections above with runnable R chunks:

```
---
title: "Incremental mediated elasticity: the P_med^delta curve"
format: html
vignette: >
  %\VignetteIndexEntry{Incremental mediated elasticity: the P_med^delta curve}
  %\VignetteEngine{quarto::html}
  %\VignetteEncoding{UTF-8}
---

```{r}
#| label: setup
#| include: false
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
set.seed(2026)
if (!requireNamespace("probmed", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
}
```
```

Use a DGP with an `A*M` interaction in the outcome so the curve is non-flat, e.g.:

```r
library(probmed)
n <- 1500
C <- rnorm(n)
A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
M <- 0.6 * A + 0.4 * C + rnorm(n)
Y <- 0.5 * A + 0.7 * M + 0.5 * A * M + 0.3 * C + rnorm(n)   # A*M interaction
d <- data.frame(A, M, Y, C)
fit <- incr_pmed(d, deltas = c(1/3, 1/2, 1, 2, 3))
fit
```

Then a sensitivity chunk:

```r
incr_sensitivity(fit, threshold = 0)
```

Write 1–3 short prose paragraphs per section (the controller's outline above is the spec for what to explain). Keep prose tight and accurate; do NOT claim double robustness for the incremental estimand (Kennedy: impossible — the estimand is itself a functional of `g`).

- [ ] **Step 2: Render the vignette to verify it builds**

Run: `R -q -e 'devtools::load_all("."); quarto::quarto_render("vignettes/incremental-pmed.qmd")'`
(or `quarto render vignettes/incremental-pmed.qmd` from the repo root)
Expected: "Output created: incremental-pmed.html", no chunk errors, the curve printout and sensitivity table appear in the output.

- [ ] **Step 3: Commit**

```bash
git add vignettes/incremental-pmed.qmd
git commit -m "docs(vignette): incremental mediated elasticity P_med^delta curve (issue #5)"
```

(Do NOT commit the rendered `.html` or `.quarto/` freeze artifacts unless the repo already tracks them for other vignettes — check `git status` and follow the existing convention.)

---

## Self-Review

**Spec coverage (Wave 2 of the design doc):** δ-sweep + reading the curve → §1–2; g-cross-fit one-step + pointwise CIs incl. term II → §3; sensitivity callout reusing Wave 1 → §4. pmed-modern reference/checkbox → out of scope (deferred, gated).

**Placeholder scan:** front-matter + guard + DGP + both calls are concrete; prose is bounded by the per-section outline.

**Consistency:** `incr_pmed()` / `incr_sensitivity()` signatures match Wave 1 and `R/incremental-pmed.R`; front-matter matches `gauge-residual.qmd` verbatim in structure.
