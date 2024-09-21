#!/usr/bin/env python3
#
# Takes a GN geno file and writes two files. One JSON file with metadata and a BIMBAM file to stdout.
#
# Note the simple JSON file contains the sample names and differs from the JSON files of GN and R/qtl2.
#
# Pjotr Prins (c) 2024

import sys
import argparse
import json
# import lmdb
# import pandas as pd
from struct import *

# import pandas as pd

meta = { "type": "gn-geno-to-gemma"
         }

parser = argparse.ArgumentParser(description='Turn GN geno format into GEMMA BIMBAM format + JSON')
parser.add_argument('file',help="GN geno file")
args = parser.parse_args()

genofn = args.file

meta["genofile"] = genofn
translate = { "A": "0", "D": "0", "B": "2", "H": "1", "U": "NA" }
header = []

with open(genofn) as f:
    for line in f:
        first = line[0]
        if first in ["#","@"]:
            header.append(line.rstrip())
            continue
        fields = line.rstrip().split("\t")
        if fields[0] == "Chr":
            # print(fields)
            samples = fields[4:]
            meta["samples"] = samples
            meta["numsamples"] = len(samples)
            continue
        # continue processing genotypes
        chrom,marker,cm,mb,*gns = fields
        # print(gns)
        outfields = [marker,"X","Y"] + [translate[item] for item in gns]
        print(",".join(outfields))

meta["header"] = header
print(json.dumps(meta),file=sys.stderr)

with open(genofn+".json","w") as f:
    f.write(json.dumps(meta))
