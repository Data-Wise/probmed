# Mediationverse CRAN Cascade — Task Checklist (2026-06-25)

> Ecosystem-wide release coordination. **Two independent arms**, joining only at the
> mediationverse meta-package: **Arm A** (medrobust → medsim, ready now) and **Arm B**
> (medfit 0.3.1 → probmed, cadence-locked to ~07-18). probmed's `incr_sensitivity`
> (PR #18, on dev) is additive — it rides probmed's pending 0.4.0 bump, not a new cascade.
> **Highest-leverage action right now: submit medrobust 0.4.0.** Authoritative state: each
> package's own `.STATUS`. Mirror in `~/mediation-planning/`.

## Dependency topology / two arms
```
Arm A (ready):   medrobust ──► medsim ──┐
Arm B (locked):  medfit 0.3.1 ──► probmed ──┤──► mediationverse (Phase 4, last)
done:            RMediation ✅ on CRAN ──────┘
off-cascade:     missingmed ✅ (r-universe, general toolkit)
```

## Version / CRAN snapshot
| Package | Local | CRAN | Gate |
|---------|-------|------|------|
| medfit | 0.3.1 (dev) | 0.2.1 ✅ | submit 0.3.1 cadence-held → ~2026-07-18 |
| RMediation | 1.5.0 | **on CRAN** ✅ | accepted/published 2026-06-20 (pkg_id 345448) — DONE |
| medrobust | 0.4.0 | — | **✅ READY TO SUBMIT** (P0): released dev→main PR #21, 0E/0W/1N. Maintainer `submit_cran()`. **Cascade bottleneck.** |
| probmed | 0.3.0.9000 | — | blocked: Imports medfit≥0.3.0 + Remotes pin |
| missingmed | 0.2.0 (released) | r-universe | **OFF the cascade** — general toolkit, not CRAN-bound now |
| medsim | 0.3.1 | — | blocked on **medrobust acceptance only**, then `submit_cran()` |
| mediationverse | 0.0.0.9000 | — | meta, Phase 4 last (needs medrobust/medsim/probmed on CRAN) |

> **Reconciled 2026-06-25** against `MediationVerse_Dashboard.md` + raw `.STATUS`. Corrected 4 stale rows
> (medrobust, RMediation, missingmed, medsim) from this doc's first draft, which had leaned on medfit's
> secondary view. Each package's **own** `.STATUS` is authoritative; cross-package mentions go stale.

---

## Phase 0 — NOW → ~2026-07-18 (parallel; cadence clock running)

### 🎯 medrobust 0.4.0 — SUBMIT NOW (highest-leverage, unblocked, P0 bottleneck)
- [ ] Build tarball from `main` (submission worktree `~/.git-worktrees/medrobust-cran`)
- [ ] `devtools::submit_cran()` (interactive + click the CRAN confirmation email)
- [ ] Flip its `.STATUS` next: → "SUBMITTED — awaiting acceptance"
- [ ] Do NOT touch `feature/cran-prep-0.4.0` worktree (superseded — reintroduces the 482s NOTE)
- Already done: revisions (PR #21, 0E/0W/1N strict `--run-donttest`), dev/main synced at 0.4.0.

### medfit (gating prerequisite for the probmed arm — but cadence-locked)
- [ ] PR `feature/fit-weights` (carries 0.3.1) → dev → main
- [ ] Finalize `cran-comments.md` for 0.3.1
- [ ] `R CMD check --as-cran` 0/0/0 — then WAIT until ~07-18 to submit (1-mo CRAN policy)

### probmed (pre-stage, can't submit until medfit 0.3.x on CRAN)
- [ ] Draft `cran-comments.md` (new submission, 5 estimators incl. incr_sensitivity)
- [ ] Confirm `--as-cran` clean **with pin still in** (baseline)
- [ ] Decide 0.4.0 NEWS content (gauge, incr+g-score, sobol, sensitivity, incr_sensitivity, wasserstein)
- [ ] Do NOT drop Remotes or bump version yet (Phase 2)

### RMediation — ✅ DONE (on CRAN, no action)
### missingmed — ✅ off the cascade (v0.2.0 released; general-toolkit roadmap is separate)

_(RMediation watcher retired — it's already accepted.)_

> **Two independent arms** (they do NOT serialize through each other):
> - **Arm A (ready now):** medrobust → medsim → mediationverse
> - **Arm B (cadence-locked):** medfit 0.3.1 → probmed
> mediationverse (Phase 4) is the single join point that waits on both arms.

---

## Arm A — medrobust → medsim (starts NOW)
**A1. medrobust** → submit_cran (Phase 0 above). Await acceptance.
**A2. medsim** (99%, gated on medrobust acceptance only — NOT probmed):
- [ ] After medrobust accepted: bump DESCRIPTION → 0.4.0 at next dev→main
- [ ] `devtools::submit_cran()`

## Arm B — medfit → probmed (cadence-locked to ~07-18)
**B1. medfit 0.3.1 → CRAN** ← critical-path event
- [ ] Submit on/after ~2026-07-18 (CRAN 1-mo policy after 0.2.1 accepted 06-18); await acceptance
**B2. probmed → 0.4.0** (after medfit 0.3.x accepted):
- [ ] Drop `Remotes: data-wise/medfit@v0.3.0`; bump 0.3.0.9000 → **0.4.0** (atomic)
- [ ] `install.packages("medfit")` (CRAN copy) → final `R CMD check --as-cran` 0/0/0
- [ ] dev → main release PR → `submit_cran()`

## Phase 4 (join) — mediationverse (meta, LAST)
- [ ] Precondition: medrobust + medsim + probmed all on CRAN (RMediation already ✅)
- [ ] Drop Remotes entries
- [ ] **Expand the thin test suite (currently 1 `test_that`)** ← dashboard flag
- [ ] `--as-cran` → submit

---

## Blockers (ordered)
1. **medrobust not yet submitted** — it's READY (PR #21); maintainer-manual `submit_cran()`. **The one unblocked bottleneck.** Do first.
2. **medfit 0.3.1 not on dev/main** — on `feature/fit-weights`, no PR. Prep now; submission still waits on #3.
3. **CRAN cadence ~2026-07-18** — hard, policy-driven, non-compressible (gates medfit→probmed arm).
4. **probmed Remotes pin** — mechanically blocked until #2 + #3 clear.

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
