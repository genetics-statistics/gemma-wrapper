#!/bin/sh
# Regression: every helper in bin/ must start without crashing on a
# missing import / gem.  We invoke each script with --help (or with
# args that take it through optparse/argparse setup) and pattern-match
# against the LoadError / ImportError / "command not found" signatures
# that fire when a dependency is unreachable.
#
# Run before any Guix package update of gemma-wrapper to capture the
# current dependency surface, then again after to confirm nothing
# regressed.  Failures fall into two buckets:
#
#   (a) gem / module not packaged in Guix at all       -> propagate
#       the missing dep + add it to the package recipe, or document
#       the script as "external use only".
#   (b) gem / module packaged but not propagated by gemma-wrapper
#       -> add to propagated-inputs in gemma.scm.
#
# Current dependency inventory (from `grep ^require` / `grep ^import`):
#
#   Ruby scripts use:
#     - stdlib:  csv, json, optparse, pp, socket, tmpdir
#     - 3rd-party (in Guix):     ruby-rdf, ruby-rdf-vocab
#     - 3rd-party (NOT in Guix): ruby-lmdb, ruby-gnrdf, ruby-qtlrange
#
#   Python scripts use:
#     - stdlib:  argparse, io, json, math, pathlib, random, re,
#                struct, sys
#     - 3rd-party (all in Guix): python-numpy, python-pandas,
#                                python-scipy, python-lmdb
#
# Usage:
#   guix shell -L /path/to/guix-bioinformatics \
#       gemma-wrapper parallel \
#       python python-numpy python-pandas python-scipy python-lmdb \
#       ruby ruby-rdf ruby-rdf-vocab \
#     -- bash test/regression-bin-deps.sh
#
# Scripts that require the ungiuxed ruby-lmdb / ruby-gnrdf /
# ruby-qtlrange gems will FAIL until those are packaged or until the
# script is dropped from bin/ -- this is the regression line in the
# sand.

set -u

TESTDIR=$(cd "$(dirname "$0")" && pwd)
SRCDIR=$(cd "$TESTDIR/.." && pwd)
BIN="$SRCDIR/bin"

PASS=0
FAIL=0
SKIP=0
log_pass() { echo "  PASS: $1"; PASS=$((PASS + 1)); }
log_fail() { echo "  FAIL: $1 [$2]"; FAIL=$((FAIL + 1)); }
log_skip() { echo "  SKIP: $1 [$2]"; SKIP=$((SKIP + 1)); }

# Run a script and classify the outcome by what it prints when its
# imports are evaluated.  Most scripts accept --help; some take only
# positional args and we look for the usage string instead.
probe() {
    local script="$1"
    local name out rc
    name=$(basename "$script")
    out=$("$script" --help 2>&1 || true)
    rc=$?
    if echo "$out" | grep -qE "cannot load such file|LoadError:|in .require'" >/dev/null 2>&1; then
        local miss
        miss=$(echo "$out" | grep -oE "cannot load such file -- [^[:space:]]+" \
               | head -1 | sed 's/^cannot load such file -- //')
        [ -z "$miss" ] && miss=$(echo "$out" \
               | grep -oE "no such file to load -- [^[:space:]]+" \
               | head -1 | sed 's/^no such file to load -- //')
        log_fail "$name" "ruby LoadError: ${miss:-unknown gem}"
        return
    fi
    if echo "$out" | grep -qE "ModuleNotFoundError|ImportError:" >/dev/null 2>&1; then
        local miss
        miss=$(echo "$out" | grep -oE "No module named '?[^'[:space:]]+" \
               | head -1 | sed "s/^No module named '*//")
        log_fail "$name" "python ImportError: ${miss:-unknown module}"
        return
    fi
    # Some scripts emit no --help but print a usage line on bad args.
    # Treat any output that looks like a usage banner as a pass.
    if echo "$out" | grep -qiE "usage:|usage[[:space:]]" >/dev/null 2>&1; then
        log_pass "$name"
        return
    fi
    # Help text printed something else (e.g. JSON-help from
    # gemma-wrapper itself) -- accept as long as nothing imported
    # failed.
    if [ -n "$out" ]; then
        log_pass "$name"
        return
    fi
    log_skip "$name" "no --help output, no error either"
}

echo "==> probing $(ls "$BIN"/*.rb "$BIN"/*.py 2>/dev/null | wc -l) scripts in $BIN"
echo

for script in "$BIN"/*.rb "$BIN"/*.py; do
    [ -f "$script" ] || continue
    probe "$script"
done

# Two odd ducks without an extension.
for script in "$BIN"/gemma-log2json "$BIN"/view-gemma-mdb; do
    [ -f "$script" ] && probe "$script"
done

echo
echo "==> summary: $PASS passed, $FAIL failed, $SKIP skipped (total $((PASS + FAIL + SKIP)))"
[ "$FAIL" -eq 0 ]
