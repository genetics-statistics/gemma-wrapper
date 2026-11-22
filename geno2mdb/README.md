# geno2mdb -- BIMBAM to pangemma's lmdb store, in Odin

The Odin counterpart of `bin/geno2mdb.py`: it turns a GEMMA BIMBAM file into the
lmdb genotype store pangemma reads with `-g X.mdb`. The store layout is
unchanged, byte for byte, because pangemma reads it:

    the lmdb environment is one file (MDB_NOSUBDIR), three databases

    "geno"    key    chr, pos, line packed big-endian as "S>L>L>"
                     (uint16 chr, uint32 pos, uint32 1-based line number)
              value  Gf: one little-endian float32 per sample, missing = NaN
                     Gb: one byte per sample, round(g*127) clipped to 0..254,
                         missing = 255
    "marker"  key    the same key -> "name<TAB>allele1<TAB>allele0", which is
                     what pangemma prints in the assoc's allele columns
    "info"    "format" -> Gf|Gb, "numsamples"/"nummarkers" (uint64 LE),
              "meta" -> {"type": "gemma-geno", "format": ..., "key-format":
              "S>L>L>", "rec-format": "f*"|"C*"}

The chromosome is a number, or the ASCII code of the letter for X/Y/M
(pangemma's CHR_X is 'X', not 23), and it is part of the key because that is how
LOCO walks the store. Keys are therefore written in chromosome and position
order, which is what the sort is for (`--no-sort` exists for inspecting a file,
not for producing a usable store).

## Why a second implementation

The python one needs an interpreter plus `numpy` and `python-lmdb` at run time.
In this deployment that meant exporting `PYTHONPATH` by hand on the API host
(otherwise: `ModuleNotFoundError: No module named 'numpy'`) and carrying three
python packages in the ravanan job's manifest for a single converter. This is
one binary and liblmdb -- and it is faster. Measured on the deployment:

| input | this binary | python | |
|---|---|---|---|
| BXD.8, 21,056 markers x 235 samples (16.5 MB) | 0.66 s | 1.21 s | 1.8x |
| HSNIH-Palmer_true, 134,918 x 6,147 (3.9 GB) | 2 min 5 s | 2 min 47 s | 1.34x |

Both stores were compared for content: identical record counts, identical
first/last keys, identical values (every record on BXD; a 1-in-97 sample plus
the counts on the 3.9 GB one). The files themselves differ by 20 KB of trailing
page padding, which is lmdb's page allocation, not the data.

## Building

Needs the Odin compiler and liblmdb. The workspace builds it with the ai-memory
development profile (which carries `odin` and `clang`) and the deployment's
`lmdb`:

    odin build src -out:geno2mdb -vet -strict-style \
      -extra-linker-flags:"-L<lmdb>/lib -Wl,-rpath,<lmdb>/lib"

`guix.scm` for this package still has to be written; it belongs next to
gemma-wrapper in guix-bioinformatics, and needs the Odin package definition
(which is not in Guix, see ai-memory's guix.scm for the shape).

## Known gaps

- **gzipped BIMBAM** is refused with a message rather than mis-parsed; linking
  zlib (`gzopen`/`gzread`) would close it, and the `.geno` input path below will
  not need it at all.
- `--geno-json` (the python's optional metadata blob) is not implemented.
- Everything is held in memory before the sort, ~3.3 GB for the panel above --
  the same as the python's `list(read_bimbam(...))`. A streaming merge-sort
  would remove that.

## Next

Reading `.geno` directly is the interesting case: it is plain text (so no zlib),
it is the file the deployment actually keeps, and the store can carry the two
things the read path currently recomputes on every request at 11.7 GB of reads
-- the sample names and a per-sample "has any genotype" mask -- in `info`.
