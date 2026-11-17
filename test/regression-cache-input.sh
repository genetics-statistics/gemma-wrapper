#!/bin/sh
# Regression tests for the genotype switches and the K that comes back out of
# the cache directory.
#
# Three things are pinned here, all of them measured on a real run:
#
#  1. A non-LOCO GWA with --input K.json and --cache-dir.  The K run archives
#     the .cXX.txt into <hash>-gemma-cXX.tar.xz and removes it from the cache
#     directory, and only the LOCO branch unpacked that archive, so GEMMA was
#     handed a -k path that did not exist and died with
#     "Enforce failed for file_kin ... open file".
#  2. The same run again must hit the cache.
#  3. LOCO with PLINK (-bfile) is refused, because GEMMA 0.98.6 ignores -loco
#     for PLINK input without saying so -- `-loco 1' and `-loco 2' produce the
#     same K -- and a wrapper that ran it would report 20 identical scans as a
#     LOCO result.
#
# Usage: test/regression-cache-input.sh [gemma-wrapper] [gemma]
#   gemma defaults to $GEMMA_COMMAND, then to `gemma' on PATH.

set -eu

TESTDIR=$(cd "$(dirname "$0")" && pwd)
SRCDIR=$(cd "$TESTDIR/.." && pwd)
BIN=${1:-$SRCDIR/bin/gemma-wrapper}
GEMMA=${2:-${GEMMA_COMMAND:-gemma}}

command -v "$GEMMA" >/dev/null 2>&1 || {
	echo "regression-cache-input: no gemma found ($GEMMA); pass it as the second argument" >&2
	exit 2
}
export GEMMA_COMMAND=$GEMMA

IN=$TESTDIR/data/input
WORK=${WORK:-$(mktemp -d -t gemma-wrapper-cache-input.XXXXXX)}
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

PASS=0
FAIL=0
log_pass() { echo "  PASS: $1"; PASS=$((PASS + 1)); }
log_fail() { echo "  FAIL: $1"; FAIL=$((FAIL + 1)); }

cd "$WORK"

echo "==> kinship with --cache-dir"
if "$BIN" --json --cache-dir "$WORK" -- \
	-g "$IN/BXD_geno.txt.gz" -p "$IN/BXD_pheno.txt" -a "$IN/BXD_snps.txt" \
	-gk >K.json 2>K.err; then
	log_pass "kinship run"
else
	log_fail "kinship run (see K.err)"
fi

echo "==> non-LOCO GWA from the cached K (the file is inside the archive)"
if "$BIN" --json --cache-dir "$WORK" --input K.json -- \
	-g "$IN/BXD_geno.txt.gz" -p "$IN/BXD_pheno.txt" -a "$IN/BXD_snps.txt" \
	-c "$IN/BXD_covariates2.txt" -lmm 2 -maf 0.1 >GWA.json 2>GWA.err; then
	if grep -q '"errno":0' GWA.json; then
		log_pass "GWA with --input from --cache-dir"
	else
		log_fail "GWA ran but reported an error"
	fi
else
	log_fail "GWA with --input from --cache-dir (see GWA.err)"
fi

echo "==> the same GWA again hits the cache"
if "$BIN" --json --cache-dir "$WORK" --input K.json -- \
	-g "$IN/BXD_geno.txt.gz" -p "$IN/BXD_pheno.txt" -a "$IN/BXD_snps.txt" \
	-c "$IN/BXD_covariates2.txt" -lmm 2 -maf 0.1 >GWA2.json 2>GWA2.err; then
	if grep -q 'CACHE HIT' GWA2.err; then
		log_pass "cache hit on the second run"
	else
		log_fail "second run recomputed instead of hitting the cache"
	fi
else
	log_fail "second GWA run"
fi

echo "==> LOCO with PLINK genotypes is refused"
# The refusal is raised before any genotype file is opened, so a prefix that
# does not exist is enough to reach it.
if "$BIN" --json --loco -- -bfile "$WORK/nonexistent" -p "$IN/BXD_pheno.txt" \
	-lmm 9 >PLINK.json 2>PLINK.err; then
	log_fail "LOCO with -bfile was accepted"
else
	if grep -q 'LOCO needs -g (BIMBAM or an lmdb .mdb store), not -bfile' PLINK.err; then
		log_pass "LOCO with -bfile refused with a reason"
	else
		log_fail "LOCO with -bfile failed for the wrong reason (see PLINK.err)"
	fi
fi

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
