#!/usr/bin/env python3
#
# Takes an R/qtl2 type CSV and writes the GEMMA file checking for the genotype/genometype names
# and writing an ordered list to stdout. Usage:
#
# rqtl2-pheno-to-gemma.py infile.csv
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

l = [] # generate genotype ordered sample output
for s in samples:
    l.append(csv.loc[s])

out = pd.DataFrame(l)
out.to_csv(sys.stdout,na_rep='NA',sep="\t",index=False,header=False)
