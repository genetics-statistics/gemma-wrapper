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
parser.add_argument('--inject',help="Inject genotypes from file")
parser.add_argument('--randomize', action=argparse.BooleanOptionalAction,help="Randomize data (default is inject haplotypes)")
parser.add_argument('--perc',type=float, default=None,help="Percentage to replace")
parser.add_argument('file',help="BIMBAM file")
args = parser.parse_args()

injectfn = args.inject
randomize = args.randomize
genofn = args.file
perc = args.perc
inject = []

if injectfn:
    with open(injectfn) as f:
        for line in f:
            fields = re.split(r"[, \t]+", line.rstrip())
            marker,x1,x2,*gns = fields
            # inject = [int(item) for item in gns]
            inject = gns
            # print(gns)

count = 0
with open(genofn) as f:
    for line in f:
        count += 1
        fields = re.split(r"[, \t]+", line.rstrip())
        marker,x1,x2,*gns = fields
        if (perc is None) or random.uniform(0,100) < perc:
            if randomize:
                gns = [str(round(random.uniform(0,2),2)) for item in gns]
            else:
                gns = inject
        outfields = [marker,x1,x2] + gns
        print(",".join(outfields))
        # if count>30:
        #     sys.exit(1)
