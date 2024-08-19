#!/usr/bin/env python3
#
# Takes an R/qtl2 type CSV and writes the GEMMA file checking for the genotype/genometype names
# from a JSON file and writing an ordered list to stdout.
#
# Note that GEMMA can take multi-column phenotype files and picks with the -n switch
# See also https://issues.genenetwork.org/topics/data/R-qtl2-format-notes
#
# Also two types of JSON file are supported. The one from GN and the one from gn-geno-to-gemma.py
#
# Pjotr Prins (c) 2024

import sys
import argparse
import json
import lmdb
import pandas as pd
from struct import *

import pandas as pd

parser = argparse.ArgumentParser(description='Turn R/qtl2 type pheno format into GEMMA pheno format')
parser.add_argument('--json',default="BXD.json",help="JSON file for sample names")
parser.add_argument('file',help="R/qtl pheno file")
args = parser.parse_args()

csv = pd.read_csv(args.file,index_col=0,na_values=['x','NA','-',''])
js = json.load(open(args.json))

if "type" in js:
    type = js["type"]
    if type == "gn-geno-to-gemma":
        samples = js["samples"]
else:
    samples = js["genofile"][0]["sample_list"] # list individuals

sampleset = set(samples)
csvset = set(csv.index)

if sampleset.difference(csvset):
    print(f"ERROR: sets differ {sampleset.difference(csvset)}",file=sys.stderr)
    sys.exit(2)

if csvset.symmetric_difference(sampleset):
    print(f"WARNING: sets differ, we'll ignore inputs from {args.file} {csvset.difference(sampleset)}",file=sys.stderr)

# if len(samples) != len(csv.index): # superfluous with above
#     print(f"ERROR: sizes do not match len(samples)={len(samples)} and len(csv.index)={len(csv.index)}",file=sys.stderr)
#     sys.exit(2)

l = [] # generate genotype ordered sample output - samples is the 'truth' list matching the genotypes - the sample has to exist on the input CSV
for s in samples:
    l.append(csv.loc[s])

out = pd.DataFrame(l)
out.to_csv(sys.stdout,na_rep='NA',sep="\t",header=False)

print(f"Wrote GEMMA pheno {len(l)} from {len(samples)} with genometypes (rows) and {len(out.columns)} collections (cols)!",file=sys.stderr)
