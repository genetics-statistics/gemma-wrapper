#!/usr/bin/env python3
#
# Takes an R/qtl2 type CSV and writes the GEMMA phenotype file
# checking for the genotype/genometype names from a JSON file, compare
# them with the first column in the CSV and write the overlap as an
# ordered list to stdout. Kind of a grep on sample names. It will also
# ERROR/WARN if samples are missing or not overlapping.
#
# Note that we use a transposed format - i.e. the individuals are the
# lines and the traits are the columns.
#
# If -n is passed in the column it checked for missing values and
# create a 'samples-reduced' list in the JSON that is picked up by
# gemma-wrapper.
#
# Note that GEMMA can take multi-column phenotype files and picks with the -n switch
# See also https://issues.genenetwork.org/topics/data/R-qtl2-format-notes
#
# Also two types of JSON file are supported. The one from GN and the one from gn-geno-to-gemma.py
#
# Pjotr Prins (c) 2024, 2025

import sys
import argparse
import json
import pandas as pd
from pathlib import Path
import math
from struct import *
import io

import pandas as pd

parser = argparse.ArgumentParser(description='Turn (transposed) R/qtl2 type pheno format into GEMMA pheno format')
parser.add_argument('--json',default="BXD.json",help="JSON file for sample names (default BXD.json)")
parser.add_argument('--header', action=argparse.BooleanOptionalAction,help="Use the header line or not")
parser.add_argument('-n',type=int, default=0,help="Get specific column number (1..) to check and skip missing values (default disabled)")
parser.add_argument('file',help="Transposed R/qtl2 pheno file (can also be json pairs)")
args = parser.parse_args()

header = 0
if args.header == False:
  header = None

filename = args.file
if ".json" in filename:
  phenojs = json.load(open(filename))
  filename = Path(filename).stem
  phenostr = f"ID,{filename}\n"
  for x in phenojs:
    line = f"{x},{phenojs[x]}\n"
    phenostr = phenostr+line
  csv = pd.read_csv(io.StringIO(phenostr),index_col=0,header=header,sep="[\s,]+",na_values=['x','NA','-',''])
  # print(phenojs)
  # print(phenostr)
else:
  csv = pd.read_csv(filename,index_col=0,header=header,sep="[\s,]+",na_values=['x','NA','-',''])

# print(csv)
# sys.exit(1)
col_n = args.n

js = json.load(open(args.json))
if "type" in js:
    type = js["type"]
    if type == "gn-geno-to-gemma":
        samples = js["samples"]
else:
    samples = js["genofile"][0]["sample_list"] # list individuals

sampleset = set(samples) # e.g. contains 'BXD1','BXD2' etc.
csvset = set(csv.index) # we take the first column to get at the sample names, also 'BXD1','BXD2' etc.

# ---- some checks
if sampleset.difference(csvset):
    print("We read the following from the CSV file:")
    print(sorted(csvset))
    print(f"ERROR: sets differ {sampleset.difference(csvset)}",file=sys.stderr)
    # print(f"\nERROR: this command has errors!")
    # sys.exit(2)

# ---- some checks
if csvset.symmetric_difference(sampleset):
    print(f"WARNING: sets differ, we'll ignore inputs from {args.file} {csvset.difference(sampleset)}",file=sys.stderr)

if col_n > 0: # get a specific column
    js["samples-column"] = col_n
    # Scan the column and collect sample names that are not NA in the samples_reduced list
    samples_reduced = {}
    # print(csv)
    header = csv.columns # trait names
    col1 = col_n-2
    trait = header[col1]
    for name,value in zip(csv.index, csv[trait]):
        if not math.isnan(value):
            # print(name,value)
            samples_reduced[name] = value
    js['trait'] = trait
    js['samples-reduced'] = samples_reduced
    js['numsamples-reduced'] = len(samples_reduced)

l = [] # generate genotype ordered sample output - samples is the 'truth' list matching the genotypes - the sample has to exist on the input CSV
for s in samples: # walk BXDs
  # if s in csv:
    l.append(csv.loc[s]) # loc gets the line
  # else:
  #   l.append("NA")

out = pd.DataFrame(l)
out.to_csv(sys.stdout,na_rep='NA',sep="\t",header=False)

outjs = filename + "-gemma.json"
with open(outjs,"w") as f:
    f.write(json.dumps(js))

# print(json.dumps(js).encode())

print(f"Wrote GEMMA pheno to stdout {len(l)} from {len(samples)} with genometypes (rows) and {len(out.columns)} collections (cols)!",file=sys.stderr)
print(f"Wrote {outjs}",file=sys.stderr)
