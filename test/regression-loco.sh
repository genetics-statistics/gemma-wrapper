#!/bin/sh
# Regression test for the LOCO + caching pipeline that gemma-wrapper
# exposes.  Mirrors the existing Rakefile test but is callable
# stand-alone (no rubygems / no bundler).  Run before any gemma update
# to capture the baseline, then re-run after the update to confirm
# nothing in the LOCO + caching path drifted.
#
# Usage:
#   guix shell -L /path/to/guix-bioinformatics \
#       gemma-wrapper parallel \
#     -- bash test/regression-loco.sh
#
# parallel must be on PATH explicitly: the gemma-wrapper Guix package
# lists it as an input but does not propagate it, so a profile that
# pulls only gemma-wrapper loses parallel at run time.  See the
# packaging TODO note at the top of test/regression-loco.sh and in
# this directory's RELEASE_NOTES.md for the fix path.
#
# Expected SHA1 hashes (from upstream README.md):
#   K hash:   1b700de28f242d561fc6769a07d88403764a996f
#   GWA hash: 9e411810ad341de6456ce0c6efd4f973356d0bad
#
# The hashes are derived from the BXD test fixtures + the GEMMA
# release pinned in the current gemma-gn2 Guix package.  When GEMMA is
# bumped, expect these to change at most in well-understood ways
# (e.g. a numerics fix); update the hashes here in lock-step with the
# package bump.

set -eu

TESTDIR=$(cd "$(dirname "$0")" && pwd)
SRCDIR=$(cd "$TESTDIR/.." && pwd)
cd "$TESTDIR"

WORK=${WORK:-$(mktemp -d -t gemma-wrapper-test.XXXXXX)}
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

K_HASH=1b700de28f242d561fc6769a07d88403764a996f
GWA_HASH=9e411810ad341de6456ce0c6efd4f973356d0bad

PASS=0
FAIL=0
log_pass() { echo "  PASS: $1"; PASS=$((PASS + 1)); }
log_fail() { echo "  FAIL: $1"; FAIL=$((FAIL + 1)); }

# Sanity: tooling on PATH.
echo "==> tool versions"
command -v gemma-wrapper >/dev/null \
    || { echo "gemma-wrapper missing on PATH"; exit 2; }
echo "    gemma-wrapper at $(command -v gemma-wrapper)"
gemma-wrapper --help 2>&1 | head -1 || true

# Bundled BXD fixtures.
GENO="$SRCDIR/test/data/input/BXD_geno.txt.gz"
PHENO="$SRCDIR/test/data/input/BXD_pheno.txt"
SNPS="$SRCDIR/test/data/input/BXD_snps.txt"
COV="$SRCDIR/test/data/input/BXD_covariates2.txt"
for f in "$GENO" "$PHENO" "$SNPS" "$COV"; do
    [ -f "$f" ] || { echo "missing fixture: $f"; exit 2; }
done

# 1. Non-LOCO -gk: compute kinship.
echo
echo "==> 1. non-LOCO kinship (-gk)"
gemma-wrapper --json --force -- \
    -g "$GENO" -p "$PHENO" -a "$SNPS" -gk -debug \
    > K0.json
grep -q "$K_HASH" K0.json \
    && log_pass "K0.json carries expected K hash" \
    || log_fail "K0.json missing K hash $K_HASH"
grep -q '"errno":0' K0.json \
    && log_pass "K0.json errno=0" \
    || log_fail "K0.json errno!=0"

# 2. Non-LOCO GWA + cache hit.
echo
echo "==> 2. non-LOCO GWA against K0.json (-lmm 2)"
gemma-wrapper --json --input K0.json -- \
    -g "$GENO" -p "$PHENO" -c "$COV" -a "$SNPS" \
    -lmm 2 -maf 0.1 -debug \
    > GWA0.json
grep -q "$GWA_HASH" GWA0.json \
    && log_pass "GWA0.json carries expected GWA hash" \
    || log_fail "GWA0.json missing GWA hash $GWA_HASH"
grep -q '"errno":0' GWA0.json \
    && log_pass "GWA0.json errno=0" \
    || log_fail "GWA0.json errno!=0"
# cache_hit on GWA0 depends on whether /tmp already cached a previous
# run, so we don't assert it here; the LOCO rerun further down does
# verify caching with a controlled state.

# 3. LOCO kinship over chromosomes 1-4.
echo
echo "==> 3. LOCO kinship (--loco --chromosomes 1,2,3,4)"
gemma-wrapper --debug --json --force \
    --loco --chromosomes 1,2,3,4 -- \
    -g "$GENO" -p "$PHENO" -a "$SNPS" -gk -debug \
    > KLOCO1.json
grep -q "$K_HASH" KLOCO1.json \
    && log_pass "KLOCO1.json carries expected K hash" \
    || log_fail "KLOCO1.json missing K hash"
grep -q '"errno":0' KLOCO1.json \
    && log_pass "KLOCO1.json errno=0" \
    || log_fail "KLOCO1.json errno!=0"

# 4. LOCO kinship again -> cache hit.
echo
echo "==> 4. LOCO kinship rerun (expect cache hit)"
gemma-wrapper --json \
    --loco --chromosomes 1,2,3,4 -- \
    -g "$GENO" -p "$PHENO" -a "$SNPS" -gk -debug \
    > KLOCO2.json
grep -q "$K_HASH" KLOCO2.json \
    && log_pass "KLOCO2.json carries expected K hash" \
    || log_fail "KLOCO2.json missing K hash"
grep -q '"cache_hit":true' KLOCO2.json \
    && log_pass "KLOCO2.json cache_hit=true" \
    || log_fail "KLOCO2.json cache_hit!=true"

# 5. LOCO GWA using cached K.
echo
echo "==> 5. LOCO GWA against KLOCO1.json"
gemma-wrapper --json --force --loco --input KLOCO1.json -- \
    -g "$GENO" -p "$PHENO" -c "$COV" -a "$SNPS" \
    -lmm 2 -maf 0.1 -debug \
    > GWA1.json
grep -q "$GWA_HASH" GWA1.json \
    && log_pass "GWA1.json carries expected GWA hash" \
    || log_fail "GWA1.json missing GWA hash"

# 6. LOCO GWA rerun -> cache hit.
echo
echo "==> 6. LOCO GWA rerun (expect cache hit)"
gemma-wrapper --json --loco --input KLOCO2.json -- \
    -g "$GENO" -p "$PHENO" -c "$COV" -a "$SNPS" \
    -lmm 2 -maf 0.1 -debug \
    > GWA2.json
grep -q "$GWA_HASH" GWA2.json \
    && log_pass "GWA2.json carries expected GWA hash" \
    || log_fail "GWA2.json missing GWA hash"
grep -q '"cache_hit":true' GWA2.json \
    && log_pass "GWA2.json cache_hit=true" \
    || log_fail "GWA2.json cache_hit!=true"

echo
echo "==> summary: $PASS passed, $FAIL failed (work dir: $WORK)"
[ "$FAIL" -eq 0 ]
