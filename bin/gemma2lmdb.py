#!/usr/bin/env python3
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

parser = argparse.ArgumentParser(description='Turn GEMMA assoc output into an lmdb db.')
parser.add_argument('--db',default="gemma.mdb",help="DB name")
parser.add_argument('--meta',required=False,help="JSON meta file name")
parser.add_argument('files',nargs='*',help="GEMMA file(s)")
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

with lmdb.open(args.db,subdir=False) as env:
    for fn in args.files:
        print(f"Processing {fn}...")
        if "log" in fn:
            with open(fn) as f:
                log[fn] = f.read()
        else:
            with open(fn) as f:
                with env.begin(write=True) as txn:
                    for line in f.readlines():
                        chr,rs,pos,miss,a1,a0,af,beta,se,l_mle,p_lrt = line.rstrip('\n').split('\t')
                        if chr=='chr':
                            continue
                        if (chr =='X'):
                            chr = X
                        elif (chr =='X'):
                            chr = Y
                        else:
                            chr = int(chr)
                        # print(f"chr={chr}, type={type(chr)}")
                        chr_c = pack('c',bytes([chr]))
                        # print(chr,chr_c,rs)
                        # key = (chr+'_'+pos).encode()
                        key = pack('>cL',chr_c,int(pos))
                        test_chr_c,test_pos = unpack('>cL',key)
                        assert chr_c == test_chr_c
                        assert test_pos == int(pos)
                        test_chr = unpack('c',chr_c)
                        # assert test_chr == int(chr), f"{test_chr} vs {int(chr)} - {chr}"
                        val = pack('=fffff', float(af), float(beta), float(se), float(l_mle), float(p_lrt))
                        assert len(val)==20, f"Packed size is expected to be 20, but is {len(val)}"
                        res = txn.put(key, bytes(val), dupdata=False, overwrite=False)
                        if res == 0:
                            print(f"WARNING: failed to update lmdb record with key {key} -- probably a duplicate {chr}:{pos} ({test_chr_c}:{test_pos})")
                        else:
                            if float(p_lrt) > 2.0:
                                hits.append([chr,int(pos),rs,p_lrt])

    with env.begin() as txn:
        with txn.cursor() as curs:
            # quick check and output of keys
            for key in list(txn.cursor().iternext(values=False)):
                chr,pos = unpack('>cL',key)
                # print(str(chr),pos)

    meta["hits"] = hits
    meta["log"] = log
    print("HELLO: ",file=sys.stderr)
    print(meta,file=sys.stderr)
    # --- Store the metadata as a JSON record in the DB
    with env.begin(write=True) as txn:
        res = txn.put('meta'.encode(), json.dumps(meta).encode(), dupdata=False, overwrite=False)
