#!/usr/bin/env python3
#
# Convert a BIMBAM mean-genotype file into the LMDB genotype store that
# pangemma reads with -g <file>.mdb.
#
# This is the Python counterpart of bin/geno2mdb.rb.  It exists because the
# ruby lmdb gem is not always installable (its C extension has to match the
# interpreter, and a gem built for one ruby is not loadable by another), while
# `python-lmdb` is packaged by Guix and Debian alike.  The on-disk format is
# the same one geno2mdb.rb writes and pangemma's io (mdb_calc_kin,
# mdb_calc_gwa) reads:
#
#   the lmdb environment is a single file (MDB_NOSUBDIR)
#
#   "geno"    key    chr, pos, line packed big-endian as "S>L>L>"
#                    (uint16 chr, uint32 pos, uint32 1-based line number)
#             value  Gf: one little-endian float32 per sample, missing = NaN
#                    Gb: one byte per sample, g*127 -> 0..254, missing = 255
#   "marker"  key    -> marker name (bytes, no allele letters)
#   "info"    "format" -> "Gb" | "Gf", "numsamples", "nummarkers", "meta"
#
# The chromosome is stored as a number, or as the ASCII code of the letter for
# X/Y/M (pangemma's CHR_X is 'X', not 23), and it is part of the key because
# that is how pangemma does leave-one-chromosome-out: it walks the "geno"
# database with a cursor and skips or stops on the chromosome in the key.  The
# keys are therefore written in chromosome and position order, not in input
# order.
#
# Note what the store does *not* hold: allele labels.  pangemma's assoc output
# reads allele1/allele0 from the marker metadata, so they come out as "-" when
# the genotypes were read from an mdb.  The genotypes themselves are the
# dosage of the first allele in the BIMBAM line (0, 1, 2), which is what
# GEMMA's BIMBAM reader uses internally as well, so results are identical.
#
# Examples:
#
#   ./bin/geno2mdb.py BXD_geno.txt.gz --snps BXD_snps.txt
#   ./bin/geno2mdb.py BXD_geno.txt.gz --snps BXD_snps.txt --format Gb
#   ./bin/geno2mdb.py --help
#
# and then, with a pangemma binary:
#
#   gemma -g BXD_geno.txt.gz.mdb -p BXD_pheno.txt -gk -o K
#   gemma -g BXD_geno.txt.gz.mdb -p BXD_pheno.txt -loco 1 -gk -o K1
#   gemma -g BXD_geno.txt.gz.mdb -p BXD_pheno.txt -loco 1 -k K1.cXX.txt -lmm 9
#
# Gf is exact (a float32 holds 0, 0.5 and 1 exactly) and is the default.  Gb is
# 4x smaller and what makes the mmap fast, but it is a quantisation: 0.5 is
# stored as 63 and read back as 63/127 = 0.496, which moves p-values in the
# third decimal.  Use Gf unless the store size is what hurts.
#
# Requires: python3 with lmdb and numpy.

import argparse
import gzip
import json
import struct
import sys

import numpy as np
import lmdb

# pangemma compares the chromosome as a uint16 in the key and prints it as a
# number, so X/Y/M are stored as their ASCII code ('X' == 88), which also sorts
# after the autosomes.
CHR_LETTERS = {"X": ord("X"), "Y": ord("Y"), "M": ord("M"), "MT": ord("M")}
KEY_FMT = ">HII"          # uint16 chr, uint32 pos, uint32 line
CHRPOS_PACK = "S>L>L>"    # the same, as named in geno2mdb.rb's metadata
MISSING_BYTE = 255


def chr_code(name):
    """The chromosome as the byte pangemma keys on, or None if unknown."""
    if name.isdigit():
        return int(name)
    return CHR_LETTERS.get(name.upper())


def open_maybe_gzip(fn):
    return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)


def read_snps(fn):
    """rs -> (chr, pos) from a GEMMA annotation file (rs pos chr per line)."""
    snps = {}
    unknown = 0
    with open(fn) as f:
        for lineno, line in enumerate(f, start=1):
            # Comma or whitespace: qtl2-impute-bimbam.R writes
            # "marker,pos,chr" while BXD_snps.txt is whitespace-separated.
            fields = line.replace(',', ' ').split()
            if len(fields) < 3:
                continue
            try:
                # R writes large positions as "9e+06", so parse a float and
                # round rather than demanding an integer.
                pos = int(round(float(fields[1])))
            except ValueError:
                raise SystemExit(
                    f"{fn}:{lineno}: expected 'rs pos chr', got {line.strip()!r}")
            code = chr_code(fields[2])
            if code is None:
                unknown += 1
                continue
            snps[fields[0]] = (code, pos)
    if unknown:
        print(f"WARNING: {unknown} markers in {fn} have a chromosome this store "
              f"cannot key (only numbers and X/Y/M)", file=sys.stderr)
    if not snps:
        raise SystemExit(f"{fn}: no usable markers -- expected 'rs pos chr' per line")
    return snps


def read_bimbam(fn, snps):
    """Yield (chr, pos, name, values, n_missing) in input order."""
    first = True
    nsamples = None
    for lineno, line in enumerate(open_maybe_gzip(fn), start=1):
        if not line.strip():
            continue
        fields = line.replace(",", " ").split()
        name = fields[0]
        if nsamples is None:
            nsamples = len(fields) - 3
        elif len(fields) - 3 != nsamples:
            raise SystemExit(f"{fn}:{lineno}: {len(fields) - 3} genotypes, "
                             f"expected {nsamples} (line 1)")
        if name not in snps:
            # A marker the annotation does not know cannot be placed in the
            # genome, and pangemma skips chromosome 0 in every LOCO run, so it
            # would be invisible.  Say so rather than store it silently.
            print(f"WARNING: marker {name} is not in the annotation, skipped",
                  file=sys.stderr)
            continue
        vals = np.array([np.nan if v in ("NA", "nan") else float(v)
                         for v in fields[3:]], dtype=np.float32)
        missing = int(np.isnan(vals).sum())
        yield snps[name][0], snps[name][1], name, vals, missing


def main():
    parser = argparse.ArgumentParser(
        description="Convert a BIMBAM mean-genotype file to the lmdb genotype "
                    "store pangemma reads (-g X.mdb).")
    parser.add_argument("files", nargs="+", metavar="BIMBAM",
                        help="BIMBAM mean genotype file(s), optionally .gz")
    parser.add_argument("--snps", required=True, metavar="FILE",
                        help="annotation: 'rs pos chr' per line "
                             "(BXD_snps.txt and friends)")
    parser.add_argument("--format", default="Gf", choices=["Gf", "Gb"],
                        help="Gf: float32, exact (default). "
                             "Gb: bytes, 4x smaller, quantises 0.5 to 63/127")
    parser.add_argument("--out", metavar="FILE",
                        help="output store (default: <input>.mdb)")
    parser.add_argument("--no-sort", action="store_true",
                        help="keep input order instead of sorting by "
                             "chromosome and position (LOCO then breaks)")
    parser.add_argument("--geno-json", metavar="FILE",
                        help="gn-geno-to-gemma style JSON to record as metadata")
    args = parser.parse_args()

    snps = read_snps(args.snps)
    if not snps:
        raise SystemExit(f"{args.snps}: no usable markers")

    meta = {"type": "gemma-geno", "format": args.format, "version": 1.0,
            "key-format": CHRPOS_PACK,
            "rec-format": "C*" if args.format == "Gb" else "f*"}
    if args.geno_json:
        with open(args.geno_json) as f:
            meta["geno"] = json.load(f)

    for fn in args.files:
        out = args.out or fn + ".mdb"
        rows = list(read_bimbam(fn, snps))
        if not rows:
            raise SystemExit(f"{fn}: no usable markers")
        nsamples = len(rows[0][3])
        total_missing = sum(r[4] for r in rows)
        if not args.no_sort:
            rows.sort(key=lambda r: (r[0], r[1]))

        env = lmdb.open(out, subdir=False, map_size=8 * 1024 ** 3, max_dbs=4)
        geno = env.open_db(b"geno")
        marker = env.open_db(b"marker")
        info = env.open_db(b"info")
        with env.begin(write=True) as txn:
            for line_no, (chr_code_, pos, name, vals, missing) in enumerate(rows, start=1):
                key = struct.pack(KEY_FMT, chr_code_, pos, line_no)
                if args.format == "Gb":
                    # Clip the *values* to 0..254 and only then put the 255
                    # sentinel in: clipping after the np.where used to clamp 255
                    # down to 254, so a missing call arrived as a genotype of
                    # 254/127 = 2.0 and nothing could tell it from a real one.
                    quantised = np.clip(np.rint(np.nan_to_num(vals, nan=0.0) * 127.0), 0,
                                        MISSING_BYTE - 1).astype(np.int32)
                    quantised = np.where(np.isnan(vals), MISSING_BYTE, quantised)
                    value = quantised.astype(np.uint8).tobytes()
                else:
                    value = vals.astype("<f4").tobytes()
                txn.put(key, value, db=geno)
                txn.put(key, name.encode(), db=marker)
            txn.put(b"format", args.format.encode(), db=info)
            txn.put(b"numsamples", struct.pack("<Q", nsamples), db=info)
            txn.put(b"nummarkers", struct.pack("<Q", len(rows)), db=info)
            txn.put(b"meta", json.dumps(meta).encode(), db=info)
        env.sync(True)
        env.close()

        pct = 100.0 * total_missing / (len(rows) * nsamples) if nsamples else 0.0
        print(f"{len(rows)} markers x {nsamples} samples, format {args.format}, "
              f"{total_missing} missing values ({pct:.1f}%) -> {out}")


if __name__ == "__main__":
    main()
