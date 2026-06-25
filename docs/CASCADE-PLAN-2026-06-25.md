# Mediationverse CRAN Cascade — Task Checklist (2026-06-25)

> Ecosystem-wide release coordination. Single gate: **medfit → CRAN**. probmed's
> `incr_sensitivity` (PR #18, on dev) is additive — it rides probmed's pending 0.4.0
> bump, not a new cascade. Authoritative state: each package's `.STATUS`.
> Mirror copy may also belong in `~/mediation-planning/`.

## Dependency topology (Imports)
```
medfit ──┬─→ probmed (>=0.3.0, PINNED) ──→ medsim ──→ mediationverse
         ├─→ missingmed (>=0.3.1, PINNED)        ↑
         ├─→ RMediation ──────────────────────────┤
         └─→ medrobust (independent) ─────────────┘
```

## Version / CRAN snapshot
| Package | Local | CRAN | Gate |
|---------|-------|------|------|
| medfit | 0.3.1 (dev) | 0.2.1 ✅ | submit 0.3.1 cadence-held → ~2026-07-18 |
| RMediation | 1.5.0 | submitted ⏳ | awaiting acceptance (pkg_id 345448) |
| medrobust | 0.4.0 | — | independent CRAN prep (6 TODOs) |
| probmed | 0.3.0.9000 | — | blocked: Imports medfit≥0.3.0 + Remotes pin |
| missingmed | 0.2.0 | — | blocked: Imports medfit≥0.3.1 |
| medsim | 0.3.1 | — | downstream (Imports medfit/probmed/RMediation/medrobust) |
| mediationverse | 0.1.0 | — | meta, last |

---

## Phase 0 — NOW → ~2026-07-18 (parallel; cadence clock running)

### medfit (gating prerequisite)
- [ ] PR `feature/fit-weights` (carries 0.3.1) → dev
- [ ] dev → main release PR; tag/prep
- [ ] Finalize `cran-comments.md` for 0.3.1
- [ ] `R CMD check --as-cran` 0/0/0

### probmed (pre-stage, can't submit yet)
- [ ] Draft `cran-comments.md` (note: new submission, 5 estimators incl. incr_sensitivity)
- [ ] Confirm `--as-cran` clean **with pin still in** (baseline)
- [ ] Decide 0.4.0 NEWS section content (gauge, incr+g-score, sobol, sensitivity, incr_sensitivity, wasserstein)
- [ ] Do NOT drop Remotes or bump version yet (Phase 2)

### medrobust 0.4.0 (independent — not gated)
- [ ] Clear 6 CRAN TODOs
- [ ] Reconcile NEWS / DESCRIPTION drift
- [ ] `--as-cran` 0/0/0 → submit (independent of medfit timing)

### RMediation (passive)
- [ ] Watch daily `medfit-cran-watch` routine for acceptance

---

## Phase 1 — ~2026-07-18: medfit 0.3.1 → CRAN  ← critical-path event
- [ ] Submit medfit 0.3.1 (earliest date per CRAN 1-month policy after 0.2.1 accepted 2026-06-18)
- [ ] Await acceptance (watcher routine)

---

## Phase 2 — after medfit 0.3.1 ACCEPTED: the pinned pair

### probmed → 0.4.0
- [ ] Drop `Remotes: data-wise/medfit@v0.3.0` from DESCRIPTION
- [ ] Bump Version 0.3.0.9000 → **0.4.0**
- [ ] Grep repo for old version strings; sync
- [ ] Final `R CMD check --as-cran` 0/0/0 (real CRAN gate now)
- [ ] dev → main release PR → submit

### missingmed
- [ ] Same unblock (Imports medfit≥0.3.1); drop Remotes, bump, `--as-cran`, submit

---

## Phase 3 — medsim
- [ ] After probmed + RMediation + medrobust all on CRAN
- [ ] Verify medfit/probmed/RMediation/medrobust pins resolve from CRAN
- [ ] `--as-cran` → submit

## Phase 4 — mediationverse (meta)
- [ ] Align all dependency versions to CRAN releases
- [ ] `--as-cran` → submit

---

## Blockers (ordered)
1. **medfit 0.3.1 not on dev/main** — stuck on `feature/fit-weights`, no PR. Resolve first.
2. **CRAN cadence ~2026-07-18** — hard, policy-driven, non-compressible.
3. **probmed Remotes pin** — mechanically blocked until 1 + 2 clear.

## Watch-outs
- probmed DESCRIPTION is 0.3.0.9000 — the 0.4.0 bump + Remotes-drop must land **together** in Phase 2, or a `.9000` dev version ships to CRAN.
- probmed is non-CRAN today by design (GitHub release pins medfit@v0.3.0); local `--as-cran` + r-universe pass only because they honor the pin — they mask the real CRAN gate.
- pmed-modern P2 manuscript submission is **independent** of this CRAN cascade (don't couple them).

---

## Appendix A — probmed Phase-2 CRAN command sequence (exact)

**PRECONDITION: medfit 0.3.1 ACCEPTED on CRAN.**

**0. Precondition check** — medfit must resolve from CRAN, not the pin:
```r
install.packages("medfit")     # pull CRAN copy
packageVersion("medfit")       # must be >= 0.3.1
```

**1. Branch + DESCRIPTION edits** (atomic — land together):
```bash
git checkout dev && git pull
git checkout -b feature/cran-0.4.0
```
In `DESCRIPTION`:
- Delete the entire `Remotes:` block (`data-wise/medfit@v0.3.0`)
- Bump `Version: 0.3.0.9000` -> `Version: 0.4.0`
- Keep `Imports: medfit (>= 0.3.0)` (CRAN 0.3.1 satisfies it)

**2. Supporting files**
- `cran-comments.md` — resubmission note (0.4.0; 5 estimators; medfit now on CRAN, Remotes dropped)
- `NEWS.md` — `# probmed 0.4.0` section (gauge, incr+g-score, sobol, sensitivity, incr_sensitivity, wasserstein)
- `grep -rn "0.3.0.9000" .` -> sync stragglers

**3. Regenerate + check** (real CRAN gate now, no pin masking):
```r
devtools::document()
rcmdcheck::rcmdcheck(args = c("--as-cran"), error_on = "warning")
# Target: 0 / 0 / 0. Prior "1 NOTE" (Remotes pin + new submission) should be gone.
```

**4. Integrate + submit**
```bash
git commit -am "release: probmed 0.4.0 (drop medfit Remotes pin; CRAN-ready)"
gh pr create --base main --head feature/cran-0.4.0 --title "Release: v0.4.0"
# after merge + CI green on main:
```
```r
devtools::submit_cran()
```

**Two traps**
1. Bump version but forget to drop Remotes -> check passes (pin masks it) but CRAN rejects.
   Verify `grep -i remotes DESCRIPTION` returns nothing.
2. Re-run `--as-cran` ONLY after `install.packages` pulls the CRAN medfit, else you test the old pinned dep.
