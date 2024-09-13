#!/usr/bin/env python3
#
# Similar to the bimbam-filter, but takes a GEMMA GRM, a genotype JSON
# file and a (pheno) file where the first column contains the
# sample/genometype/strain names only the subset will be picked from
# the GRM file and written out. E.g.
#
#   ./bin/grm-filter.py --json BXD.geno.json --sample-file BXD_pheno_Dave-GEMMA-samples.txt output/result.cXX.txt
#
# Note the simple JSON file from gn-geno-to-gemma.py contains the
# sample names and differs from the JSON files of GN and R/qtl2. This
# script outputs a similar file.
#
# Pjotr Prins (c) 2024

import sys
import argparse
import json
import re

meta = { "type": "grm-filter"
         }

parser = argparse.ArgumentParser(description='Turn GRM into a filtered/reduced GEMMA GRM')
parser.add_argument('--json',required=True,help="JSON file contains names of GRM samples/genometypes")
parser.add_argument('--sample-file',required=True,help="Use sample/pheno file to filter for genometypes/samples using the 1st column")
parser.add_argument('file',help="GEMMA GRM file")
args = parser.parse_args()

# Fetch GRM samples
json_in = json.load(open(args.json))
if 'samples' in json_in:
    all_samples = json_in['samples']
else:
    # picks up gemma-wrapper K.json
    all_samples = json_in['input']['samples']
count_all_samples = len(all_samples)

# Get samples to filter
filter_samples = []
with open(args.sample_file) as f:
    for line in f.readlines():
        sample = line.rstrip().split("\t")[0]
        filter_samples.append(sample)
# print(filter_samples)

# Create an index on samples
index = []
for s in filter_samples:
    index.append(all_samples.index(s))

grmfn = args.file

num = 0
with open(grmfn) as f:
    for line in f:
        if num in index:
            fields = re.split(r"[, \t]+", line.rstrip())
            # print(fields)
            assert len(fields) == count_all_samples, f"{len(fields)} != {count_all_samples}"
            # marker,x1,x2,*gns = fields
            # print(fields)
            filtered = [fields[i] for i in index]
            # [translate[item] for item in gns
            # filtered = [gns[i] for i in index]
            # outfields = [marker,x1,x2] + filtered
            print("\t".join(filtered))
        num += 1

# meta["header"] = header
# print(json.dumps(meta),file=sys.stderr)

# with open(genofn+".json","w") as f:
#     f.write(json.dumps(meta))
