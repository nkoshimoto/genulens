#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT_DIR"

if [[ ! -x backup/pre-refactor-bin/genulens ]]; then
  echo "backup/pre-refactor-bin/genulens missing; skipping regression comparison"
  exit 0
fi

backup/pre-refactor-bin/genulens Nsimu 100 seed 1234 >/tmp/genulens-old.out
./genulens Nsimu 100 seed 1234 >/tmp/genulens-new.out

grep -E '# Nlike/N/Ngen=' /tmp/genulens-old.out >/tmp/genulens-old-summary.out
grep -E '# Nlike/N/Ngen=' /tmp/genulens-new.out >/tmp/genulens-new-summary.out
diff -u /tmp/genulens-old-summary.out /tmp/genulens-new-summary.out

echo "regression CLI tests passed"

