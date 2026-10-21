#!/bin/sh
# Exercise publication guards without running R or writing live GitHub statuses.
set -eu
gate=$(cd "$(dirname "$0")" && pwd)/check-full-local.sh
fixture=$(mktemp -d "${TMPDIR:-/tmp}/probmed-local-ci-test.XXXXXX")
trap 'rm -r "$fixture"' EXIT
mkdir "$fixture/bin" "$fixture/repo"
cp "$gate" "$fixture/repo/gate.sh"
printf 'RoxygenNote: 7.3.3\n' > "$fixture/repo/DESCRIPTION"
cat > "$fixture/bin/gh" <<'STUB'
#!/bin/sh
set -eu
case "$*" in
  *headRefOid*)
    if [ "$CASE" = moved ] && [ -e "$FIXTURE/tested" ]; then echo moved; else git rev-parse HEAD; fi ;;
  *baseRefName*) echo dev ;;
  *isCrossRepository*) if [ "$CASE" = fork ]; then echo true; else echo false; fi ;;
  *--json\ state*) echo OPEN ;;
  api*) printf '%s\n' "$*" >> "$FIXTURE/statuses" ;;
  *) exit 1 ;;
esac
STUB
cat > "$fixture/bin/Rscript" <<'STUB'
#!/bin/sh
set -eu
case "$*" in
  *devtools::test*)
    touch "$FIXTURE/tested"
    case "$CASE" in
      dirty) printf 'changed\n' >> DESCRIPTION ;;
      head) git -c user.name=Fixture -c user.email=fixture@example.invalid commit --quiet --allow-empty -m moved ;;
      stage-failure) exit 1 ;;
    esac ;;
  *pkgdown::*)
    for site do :; done
    mkdir -p "$site/docs/reference"
    printf 'index\n' > "$site/docs/index.html"
    printf 'reference\n' > "$site/docs/reference/index.html" ;;
esac
STUB
chmod +x "$fixture/bin/gh" "$fixture/bin/Rscript"
export FIXTURE="$fixture" PATH="$fixture/bin:$PATH"
cd "$fixture/repo"
git init -q
git config user.name Fixture
git config user.email fixture@example.invalid
git remote add origin https://github.com/Data-Wise/probmed.git
git add .
git commit -qm fixture
initial=$(git rev-parse HEAD)
count=0
for CASE in success dirty head moved stage-failure fork; do
  export CASE
  git switch --quiet --detach "$initial"
  git restore DESCRIPTION
  rm -f "$fixture/statuses" "$fixture/tested"
  result=0
  sh ./gate.sh --pr 53 > "$fixture/output" 2>&1 || result=$?
  if [ "$CASE" = success ]; then
    [ "$result" -eq 0 ]
    grep -q 'state=success' "$fixture/statuses"
  else
    [ "$result" -ne 0 ]
    if [ "$CASE" = fork ]; then
      [ ! -e "$fixture/statuses" ]
    else
      grep -q 'state=failure' "$fixture/statuses"
      if grep -q 'state=success' "$fixture/statuses"; then exit 1; fi
    fi
  fi
  count=$((count + 1))
done
printf 'PASS: %s local-CI publication scenarios\n' "$count"
