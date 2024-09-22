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
import os
from pathlib import Path
import re
from struct import *
import tempfile

from signal import signal, SIGPIPE, SIG_DFL # avoid broken pipe error
signal(SIGPIPE,SIG_DFL)

parser = argparse.ArgumentParser(description="Fetch GEMMA lmdb values.")
parser.add_argument('--anno',required=False,help="SNP annotation file with the format 'rs31443144, 3010274, 1'")
parser.add_argument('--sort',action=argparse.BooleanOptionalAction,default=True,help="Sort on significance")
parser.add_argument('lmdb',nargs='?',help="GEMMA lmdb db file name (also can take tar.xz)")
args = parser.parse_args()

# ASCII
X=ord('X')
Y=ord('Y')

snps = {}
if args.anno:
    for line in open(args.anno, 'r'):
        snp,pos,chrom = re.split(r"\s+",line.rstrip('\n'))
        key = chrom+":"+pos
        snps[key] = snp

print("chr,pos,marker,af,beta,se,l_mle,l_lrt,-logP")

with tempfile.TemporaryDirectory() as tmpdir:
    print(f"Created temporary directory {tmpdir}",file=sys.stderr)

    fn = args.lmdb
    if fn.endswith('.xz'):
        print(f"Unpack {fn}...",file=sys.stderr)
        os.system(f"tar xvJf {fn} -C {tmpdir} > /dev/null")
        fn = tmpdir+"/"+Path(Path(fn).stem).stem+".mdb"

    print(f"Reading {fn}...",file=sys.stderr)

    result = []
    with lmdb.open(fn,subdir=False) as env:
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
                        result.append([chr,str(pos),snp,str(round(af,4)),str(round(effect,4)),
                                       str(round(se,4)),str(round(l_mle,4)),str(round(p_lrt,4)),str(round(minusLogP,2))])

if args.sort:
    result = sorted(result, key=lambda x: x[8], reverse=True)

for l in result:
    print(",".join(l))
    #      print(",".join([chr,str(pos),snp,str(round(af,4)),str(round(effect,4)),
    #              str(round(se,4)),str(round(l_mle,4)),str(round(p_lrt,4)),str(round(minusLogP,2))]))
