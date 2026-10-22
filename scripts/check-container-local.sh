#!/bin/sh
# Linux R CMD check in dedicated containers (R release, oldrel-1, devel).
# Commit-bound: checks `git archive HEAD`, never the working tree. Use
# --pr NUMBER to publish the result as the local-ci/probmed-container status.
# Windows and macOS are not covered here; the hosted R-CMD-check workflow keeps them.
set -eu

repo=Data-Wise/probmed
context=local-ci/probmed-container
pr=
tags=${PROBMED_R_TAGS:-"4.5 4.4 devel"}
if [ "$#" -ne 0 ]; then
  [ "$#" -eq 2 ] && [ "$1" = "--pr" ] || {
    echo 'Usage: scripts/check-container-local.sh [--pr NUMBER]' >&2
    exit 2
  }
  pr=$2
  case "$pr" in ''|*[!0-9]*) echo 'PR number must be numeric' >&2; exit 2 ;; esac
fi

command -v docker >/dev/null || { echo 'docker is required.' >&2; exit 2; }
docker info >/dev/null 2>&1 || { echo 'The Docker daemon is not running.' >&2; exit 2; }

root=$(git rev-parse --show-toplevel)
cd "$root"
if [ -n "$(git status --porcelain --untracked-files=normal)" ]; then
  echo 'Container CI requires a clean checkout of the commit to be tested.' >&2
  exit 2
fi
sha=$(git rev-parse HEAD)
work=$(mktemp -d "${TMPDIR:-/tmp}/probmed-container.XXXXXX")
work=$(cd "$work" && pwd -P)
echo "Container CI for $sha; logs: $work"

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
  pr_state=$(gh pr view "$pr" --repo "$repo" --json state --jq .state)
  pr_fork=$(gh pr view "$pr" --repo "$repo" --json isCrossRepository --jq .isCrossRepository)
  [ "$pr_sha" = "$sha" ] && [ "$pr_state" = OPEN ] || {
    echo 'PR is not open at this exact local commit.' >&2; exit 2;
  }
  [ "$pr_fork" = false ] || { echo 'Refusing to run fork PR code locally.' >&2; exit 2; }
  post_status pending 'Linux container gate running'
fi

finish() {
  result=$?
  if [ "$result" -ne 0 ]; then
    post_status failure 'Linux container gate failed; inspect local logs' || true
    echo "Container CI FAILED; logs: $work" >&2
  fi
}
trap finish EXIT

mkdir "$work/src"
git archive HEAD | tar -x -C "$work/src"

for tag in $tags; do
  echo "Building image R $tag"
  docker build -q -f scripts/docker/Dockerfile --build-arg "R_TAG=$tag" \
    -t "probmed-ci:$tag" "$work/src" >"$work/build-$tag.log" 2>&1 || {
    tail -n 25 "$work/build-$tag.log" >&2; exit 1; }
  echo "Checking on R $tag"
  docker run --rm -e NOT_CRAN=true -v "$work/src:/pkg:ro" "probmed-ci:$tag" \
    Rscript -e 'x <- rcmdcheck::rcmdcheck("/pkg", args = "--no-manual", error_on = "never"); print(x); if (length(x$errors) || length(x$warnings)) quit(status = 1)' \
    >"$work/check-$tag.log" 2>&1 || {
    tail -n 40 "$work/check-$tag.log" >&2; exit 1; }
done

if [ "$(git rev-parse HEAD)" != "$sha" ] || [ -n "$(git status --porcelain --untracked-files=normal)" ]; then
  echo 'Checkout changed during container CI; refusing success.' >&2
  exit 1
fi
if [ -n "$pr" ]; then
  latest=$(gh pr view "$pr" --repo "$repo" --json headRefOid --jq .headRefOid)
  [ "$latest" = "$sha" ] || { echo 'PR head moved during container CI; refusing success.' >&2; exit 1; }
  post_status success "Linux R CMD check passed on R $tags at this commit"
fi
trap - EXIT
echo "Container CI PASSED for $sha; logs: $work"
