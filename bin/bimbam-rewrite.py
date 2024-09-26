#!/usr/bin/env python3
#
# Takes a BIMBAM file and rewrites it with modified haplotypes.
#
# Pjotr Prins (c) 2024

import sys
import argparse
import json
import random
import re
# from struct import *

# import pandas as pd

parser = argparse.ArgumentParser(description='Rewrite BIMBAM')
parser.add_argument('file',help="BIMBAM file")
args = parser.parse_args()

genofn = args.file

with open(genofn) as f:
    for line in f:
        fields = re.split(r"[, \t]+", line.rstrip())
        # print(fields)
        # assert len(fields) == count_all_samples+3, f"{len(fields)} != {count_all_samples+3}"
        marker,x1,x2,*gns = fields
        # print(gns)
        gns = [str(round(random.uniform(0,2),2)) for item in gns]
        # print(gns)
        outfields = [marker,x1,x2] + gns
        print(",".join(outfields))
        # sys.exit(1)
