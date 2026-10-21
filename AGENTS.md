# AGENTS.md for probmed Package

## Working contract

P_med is tie-adjusted stochastic superiority: P(Y(x,M(x)) > Y(x,M(x*))) + half the tie probability, with treatment held fixed and only mediator changing. Preserve cross-world assumptions and distinguish estimand definition from implemented coupling/model support. Do not conflate ties, strict benefit, or a null score with absence of every mediation mechanism. Gaussian and Bernoulli support does not establish ordinal/count support.
Extraction belongs to medfit; probmed owns its bootstrap and PmedResult. Preserve family/link-scale simulation, direct-effect cancellation, bounded scores and reproducible bootstrap. Read Code Architecture and Testing Strategy before computation/inference work; preserve the MBCO union-null branches and scientific verification requirements.
R >=4.1, S7, native pipe, explicit namespaces, snake_case/CamelCase conventions. Required: document/test, `rcmdcheck::rcmdcheck(args="--no-manual", error_on="error")`, pkgdown verification; coverage >90% overall, 100% computation core. Test oracle accuracy, c'-invariance, binary non-degeneracy, extraction and edge cases. main <- dev <- feature worktrees; retain current release/submission holds.

## On-demand requirements

The exact pre-condensation instructions are preserved in `/Users/davood/.codex/reference/repos/r-packages-active-probmed/AGENTS-detail-2026-10-12.md`. Before a task, read only its applicable sections: scientific/model changes, builds/tests, source conventions, dependency/API changes, packaging/releases, configuration/imports, deployment, or workflow cleanup. Every applicable requirement and verification gate there remains binding; this shorter file does not waive any. Do not load the whole reference for unrelated tasks. Historical versions, counts, results and completed-work snapshots are not live status: check the owning source and .STATUS. Current user instructions and project holds take precedence.

Preserve unrelated changes and concurrent sessions. Required checks must pass before claiming completion; inspect actual rendered artifacts where applicable. Author review/submission holds remain active until explicitly lifted. Passing checks never authorizes submission, merge, release, deployment, email sending or privacy changes.
