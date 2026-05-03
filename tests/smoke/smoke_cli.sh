#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT_DIR"

./genulens --help >/tmp/genulens-help.out
grep -q '^Usage: ./genulens ' /tmp/genulens-help.out

./pre_gapmoe/calc_rho_profile --help >/tmp/calc-rho-help.out
grep -q 'Usage:' /tmp/calc-rho-help.out

./pre_gapmoe/calc_mass_dist --help >/tmp/calc-mass-help.out
grep -q 'Usage:' /tmp/calc-mass-help.out

./pre_gapmoe/calc_murel_dist --help >/tmp/calc-murel-help.out
grep -q 'Usage:' /tmp/calc-murel-help.out

tmpdir="$(mktemp -d)"
"$ROOT_DIR/genulens" l 0.5 b 0.2 Nsimu 3 seed 1234 >"$tmpdir/genulens-from-other-cwd.out"
grep -q '# Nlike/N/Ngen=' "$tmpdir/genulens-from-other-cwd.out"
rm -rf "$tmpdir"

echo "smoke CLI tests passed"
