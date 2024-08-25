#!/usr/bin/env python3
#
# Takes a BIMBAM file, a JSON file and a (pheno) file where the first
# column contains the sample/genometype/strain names only the subset
# will be picked from the genotype file and written out. E.g.
#
#   grep -v "NA" BXD_pheno_Dave-GEMMA.txt > BXD_pheno_Dave-GEMMA-samples.txt
#   ./bin/bimbam-filter.py --json BXD.geno.json --sample-file BXD_pheno_Dave-GEMMA-samples.txt BXD_geno.txt

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
# import lmdb
# import pandas as pd
# from struct import *

# import pandas as pd

meta = { "type": "bimbam-filter"
         }

parser = argparse.ArgumentParser(description='Turn BIMBAM into a filtered GEMMA BIMBAM format + JSON')
parser.add_argument('--json',required=True,help="JSON file contains names of samples/genometypes")
parser.add_argument('--sample-file',required=True,help="Use sample/pheno file to filter for genometypes/samples using the 1st column")
parser.add_argument('file',help="BIMBAM file")
args = parser.parse_args()

# Fetch BIMBAM samples
json_in = json.load(open(args.json))
all_samples = json_in['samples']
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
i = 0
for s in filter_samples:
    index.append(all_samples.index(s))

genofn = args.file
# meta["genofile"] = genofn

with open(genofn) as f:
    for line in f:
        fields = re.split(r"[, \t]+", line.rstrip())
        # print(fields)
        assert len(fields) == count_all_samples+3, f"{len(fields)} != {count_all_samples+3}"
        marker,x1,x2,*gns = fields
        # print(gns)
        # [translate[item] for item in gns
        filtered = [gns[i] for i in index]
        outfields = [marker,x1,x2] + filtered
        print(",".join(outfields))

# meta["header"] = header
# print(json.dumps(meta),file=sys.stderr)

# with open(genofn+".json","w") as f:
#     f.write(json.dumps(meta))
