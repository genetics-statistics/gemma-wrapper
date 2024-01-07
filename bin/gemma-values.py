#!/usr/bin/env python3
#
# Get values from the gemma-lmdb file.
#
# Note that the chr+position is used in the key where position is
# stored big-endian to allow for proper sorting(!) X and Y chromosomes
# are stored as their ASCII value (88 and 89(.
#
# The records contain the standard gemma output stored as floats.
#
# The conversion will fail if:
#
# - two markers share the same position

import sys
import argparse
import json
import lmdb
import math
from struct import *

parser = argparse.ArgumentParser(description="Fetch GEMMA lmdb values.")
parser.add_argument('--anno',required=False,help="SNP annotation file with the format 'rs31443144, 3010274, 1'")
parser.add_argument('lmdb',nargs='?',help="GEMMA lmdb db file name")
args = parser.parse_args()

# ASCII
X=ord('X')
Y=ord('Y')

snps = {}
if args.anno:
    for line in open(args.anno, 'r'):
        snp,pos,chrom = line.rstrip('\n').split(", ")
        key = chrom+":"+pos
        snps[key] = snp

print("chr,pos,marker,af,beta,se,l_mle,-logP")
with lmdb.open(args.lmdb,subdir=False) as env:
    with env.begin() as txn:
        with txn.cursor() as curs:
            # quick check and output of keys
            for key in list(txn.cursor().iternext(values=False)):
                if key != b'meta':
                    chr1,pos = unpack('>cL',key)
                    chr2 = int.from_bytes(chr1,"little")
                    # chr2 = chr1
                    if chr2 == X:
                        chr = "X"
                    elif chr1 == Y:
                        chr = "Y"
                    else:
                        chr = str(chr2)
                    rec = txn.get(key)
                    af,beta,se,l_mle,p_lrt = unpack('=fffff',rec)
                    effect = -beta/2.0
                    minusLogP = -math.log(p_lrt,10)
                    snp = "?"
                    if args.anno:
                        key2 = chr+":"+str(pos)
                        snp = snps[key2]

                    print(",".join([chr,str(pos),snp,str(round(af,4)),str(round(effect,4)),str(round(se,4)),str(l_mle),str(round(minusLogP,2))]))
