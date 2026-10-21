#!/bin/sh
# Full, commit-bound local gate. Use --pr NUMBER to publish its result to GitHub.
# shellcheck disable=SC2016 # R expressions below deliberately use single quotes.
set -eu

repo=Data-Wise/probmed
context=local-ci/probmed-full
pr=
if [ "$#" -ne 0 ]; then
  [ "$#" -eq 2 ] && [ "$1" = "--pr" ] || {
    echo 'Usage: scripts/check-full-local.sh [--pr NUMBER]' >&2
    exit 2
  }
  pr=$2
  case "$pr" in ''|*[!0-9]*) echo 'PR number must be numeric' >&2; exit 2 ;; esac
fi

root=$(git rev-parse --show-toplevel)
cd "$root"
if [ -n "$(git status --porcelain --untracked-files=normal)" ]; then
  echo 'Full local CI requires a clean checkout of the commit to be tested.' >&2
  exit 2
fi
sha=$(git rev-parse HEAD)
work=$(mktemp -d "${TMPDIR:-/tmp}/probmed-full-local.XXXXXX")
echo "Local CI for $sha; logs: $work"

post_status() {
  [ -n "$pr" ] || return 0
  gh api --silent --method POST "repos/$repo/statuses/$sha" \
    -f "state=$1" -f "context=$context" -f "description=$2"
}
if [ -n "$pr" ]; then
  command -v gh >/dev/null || { echo 'gh CLI is required for --pr.' >&2; exit 2; }
  remote=$(git remote get-url origin)
  case "$remote" in
    https://github.com/Data-Wise/probmed.git|git@github.com:Data-Wise/probmed.git) ;;
    *) echo "Unexpected origin: $remote" >&2; exit 2 ;;
  esac
  pr_sha=$(gh pr view "$pr" --repo "$repo" --json headRefOid --jq .headRefOid)
  pr_base=$(gh pr view "$pr" --repo "$repo" --json baseRefName --jq .baseRefName)
  pr_state=$(gh pr view "$pr" --repo "$repo" --json state --jq .state)
  pr_fork=$(gh pr view "$pr" --repo "$repo" --json isCrossRepository --jq .isCrossRepository)
  [ "$pr_sha" = "$sha" ] && [ "$pr_state" = OPEN ] || {
    echo 'PR is not open at this exact local commit.' >&2; exit 2;
  }
  [ "$pr_fork" = false ] || { echo 'Refusing to run fork PR code locally.' >&2; exit 2; }
  case "$pr_base" in main|dev) ;; *) echo 'PR base must be main or dev.' >&2; exit 2 ;; esac
  post_status pending 'Full local R gate running'
fi

finish() {
  result=$?
  if [ "$result" -ne 0 ]; then
    post_status failure 'Full local R gate failed; inspect local logs' || true
    echo "Full local CI FAILED; logs: $work" >&2
  fi
}
trap finish EXIT

stage() {
  label=$1
  shift
  echo "Running $label"
  if ! "$@" >"$work/$label.log" 2>&1; then
    tail -n 25 "$work/$label.log" >&2
    return 1
  fi
}

roxygen_version=$(sed -n 's/^RoxygenNote: *//p' DESCRIPTION)
[ -n "$roxygen_version" ] || { echo 'DESCRIPTION must record RoxygenNote.' >&2; exit 1; }
stage roxygen-version Rscript -e "if (as.character(utils::packageVersion('roxygen2')) != '$roxygen_version') stop('Expected roxygen2 $roxygen_version; use an isolated R library')"
stage document Rscript -e 'devtools::document()'
if [ -n "$(git status --porcelain --untracked-files=normal)" ]; then
  echo 'Documentation changed; review and commit generated files, then rerun.' >&2
  exit 1
fi
export NOT_CRAN=true
stage tests Rscript -e 'devtools::test(stop_on_failure = TRUE)'
stage rcmdcheck Rscript -e 'x <- rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "error"); if (length(x$warnings)) stop("R CMD check warnings: ", paste(x$warnings, collapse = "; "))'
stage coverage Rscript -e 'cov <- covr::package_coverage(path = ".", type = "tests", quiet = TRUE, clean = FALSE); overall <- covr::percent_coverage(cov); tally <- covr::tally_coverage(cov, by = "expression"); core <- tally[grepl("(^|/)R/compute-core[.]R$", tally$filename), , drop = FALSE]; cat("Overall:", overall, "core expressions:", nrow(core), "\n"); if (!is.finite(overall) || overall <= 90 || !nrow(core) || anyNA(core$value) || any(core$value <= 0)) stop("Coverage targets not met")'

mkdir "$work/site-source"
stage site-copy rsync -a --exclude=.git --exclude=docs --exclude=AGENTS.md --exclude=CLAUDE.md "$root/" "$work/site-source/"
stage pkgdown Rscript -e 'pkgdown::build_site_github_pages(pkg = commandArgs(trailingOnly = TRUE)[1], new_process = FALSE, install = TRUE)' "$work/site-source"
site="$work/site-source/docs"
if [ -s "$site/dev/index.html" ]; then site="$site/dev"; fi
[ -s "$site/index.html" ] && [ -s "$site/reference/index.html" ] || {
  echo 'pkgdown index or reference page missing' >&2; exit 1;
}

if [ "$(git rev-parse HEAD)" != "$sha" ] || [ -n "$(git status --porcelain --untracked-files=normal)" ]; then
  echo 'Checkout changed during local CI; refusing success.' >&2
  exit 1
fi
if [ -n "$pr" ]; then
  latest=$(gh pr view "$pr" --repo "$repo" --json headRefOid --jq .headRefOid)
  [ "$latest" = "$sha" ] || { echo 'PR head moved during local CI; refusing success.' >&2; exit 1; }
  post_status success 'Full local R gate passed at this commit'
fi
trap - EXIT
echo "Full local CI PASSED for $sha; logs: $work"
