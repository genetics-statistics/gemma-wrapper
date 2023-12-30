#!/usr/bin/env python3

import sys
import argparse
import lmdb
from struct import *

parser = argparse.ArgumentParser(description='Turn GEMMA assoc output into an lmdb db.')
parser.add_argument('--db',default="gemma",help="DB name")
parser.add_argument('files',nargs='*',help="GEMMA assoc file(s)")
args = parser.parse_args()

with lmdb.open(args.db+".mdb",subdir=False) as env:
    for fn in args.files:
        print(f"Processing {fn}...")
        with open(fn) as f:
            with env.begin(write=True) as txn:
                for line in f.readlines():
                    chr,rs,pos,miss,a1,a0,af,logl_H1,l_mle,p_lrt = line.rstrip('\n').split('\t')
                    if chr=='chr':
                        continue
                    print(chr,rs)
                    key = (chr+'_'+pos).encode()
                    key = pack('=2sL',chr.encode(),int(pos))
                    val = pack('=ffff', float(af), float(logl_H1), float(l_mle), float(p_lrt))
                    assert len(val)==16, f"Packed size is expected to be 16, but is {len(val)}"
                    res = txn.put(key, bytes(val), dupdata=False, overwrite=False)
                    assert res,f"Failed to update lmdb record with key {key} -- probably a duplicate"
