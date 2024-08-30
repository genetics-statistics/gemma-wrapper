#!/usr/bin/env python3
#
# Takes an R/qtl2 type CSV and writes the GEMMA file checking for the
# genotype/genometype names from a JSON file, compare them with the
# first column in the CSV and write the overlap as an ordered list to
# stdout. Kind of a grep on sample names. It will also ERROR/WARN if
# samples are missing or not overlapping.
#
# If -n is passed in the column it checked for missing values and drops those samples.
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
# import lmdb
import pandas as pd
import math
from struct import *

import pandas as pd

parser = argparse.ArgumentParser(description='Turn R/qtl2 type pheno format into GEMMA pheno format')
parser.add_argument('--json',default="BXD.json",help="JSON file for sample names (default BXD.json)")
parser.add_argument('--header', action=argparse.BooleanOptionalAction,help="Use the header line or not")
parser.add_argument('-n',type=int, default=0,help="Column number to check and skip missing values (default 1)")
parser.add_argument('file',help="R/qtl pheno file")
args = parser.parse_args()

header = 0
if args.header == False:
  header = None

csv = pd.read_csv(args.file,index_col=0,header=header,sep="[\s,]+",na_values=['x','NA','-',''])
col = args.n

js = json.load(open(args.json))
if "type" in js:
    type = js["type"]
    if type == "gn-geno-to-gemma":
        samples = js["samples"]
else:
    samples = js["genofile"][0]["sample_list"] # list individuals

sampleset = set(samples)
csvset = set(csv.index) # we take the first column to get at the sample names

if sampleset.difference(csvset):
    print("We read the following from the CSV file:")
    print(sorted(csvset))
    print(f"ERROR: sets differ {sampleset.difference(csvset)}",file=sys.stderr)
    sys.exit(2)

if csvset.symmetric_difference(sampleset):
    print(f"WARNING: sets differ, we'll ignore inputs from {args.file} {csvset.difference(sampleset)}",file=sys.stderr)

if col > 0:
    js["samples-column"] = col
    # Walk the column
    samples_reduced = {}
    for name,value in zip(csv.index, csv[1]):
        if not math.isnan(value):
            print(name,value)
            samples_reduced[name] = value
    js['samples-reduced'] = samples_reduced
    js['numsamples-reduced'] = len(samples_reduced)

l = [] # generate genotype ordered sample output - samples is the 'truth' list matching the genotypes - the sample has to exist on the input CSV
for s in samples:
    l.append(csv.loc[s])

out = pd.DataFrame(l)
out.to_csv(sys.stdout,na_rep='NA',sep="\t",header=False)

outjs = args.file + ".json"
with open(outjs,"w") as f:
    f.write(json.dumps(js))

# print(json.dumps(js).encode())

print(f"Wrote GEMMA pheno {len(l)} from {len(samples)} with genometypes (rows) and {len(out.columns)} collections (cols)!",file=sys.stderr)
print(f"Wrote {outjs}")
