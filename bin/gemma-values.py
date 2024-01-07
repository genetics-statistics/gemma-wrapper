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
from struct import *

parser = argparse.ArgumentParser(description="Fetch GEMMA lmdb values.")
# parser.add_argument('--db',default="gemma.mdb",help="DB name")
# parser.add_argument('--meta',required=False,help="JSON meta file name")

parser.add_argument('lmdb',nargs='?',help="GEMMA lmdb db file name")
args = parser.parse_args()

# ASCII
X=ord('X')
Y=ord('Y')

meta = { "type": "gemma-assoc",
         "version": 1.0,
         "key-format": ">cL",
         "rec-format": "=ffff" }
log = [] # track output log
hits = [] # track hits

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
                    print(chr,pos)
