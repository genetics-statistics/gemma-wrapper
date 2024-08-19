#!/usr/bin/env python3
#
# Takes an R/qtl2 type CSV and writes the GEMMA file checking for the genotype/genometype names
# and writing an ordered list to stdout. Usage:
#
# Note that GEMMA can take multi-column phenotype files and picks with the -n switch
# See also https://issues.genenetwork.org/topics/data/R-qtl2-format-notes
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
parser.add_argument('file',help="R/qtl pheno file")
args = parser.parse_args()

csv = pd.read_csv(args.file,index_col=0,na_values=['x','NA','-',''])
js = json.load(open("BXD.json"))
samples = js["genofile"][0]["sample_list"] # list individuals

sampleset = set(samples)
csvset = set(csv.index)

if sampleset.difference(csvset):
    print(f"ERROR: sets differ {sampleset.difference(csvset)}",file=sys.stderr)
    sys.exit(2)

if csvset.difference(sampleset):
    print(f"WARNING: sets differ, we'll ignore inputs from {args.file} {csvset.difference(sampleset)}",file=sys.stderr)

# if len(samples) != len(csv.index): # superfluous with above
#     print(f"ERROR: sizes do not match len(samples)={len(samples)} and len(csv.index)={len(csv.index)}",file=sys.stderr)
#     sys.exit(2)

l = [] # generate genotype ordered sample output
for s in samples:
    l.append(csv.loc[s])

out = pd.DataFrame(l)
out.to_csv(sys.stdout,na_rep='NA',sep="\t",index=False,header=False)

print(f"Wrote GEMMA pheno {len(l)} rows and {len(out.columns)} cols!",file=sys.stderr)
