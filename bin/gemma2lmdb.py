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
import numpy as np
import scipy.stats as stat
import math

SIGNIFICANT = 4.0

parser = argparse.ArgumentParser(description='Turn GEMMA assoc output into an lmdb db.')
parser.add_argument('--db',default="gemma.mdb",help="DB name")
parser.add_argument('--meta',required=False,help="JSON meta file name")
parser.add_argument('files',nargs='*',help="GEMMA file(s)")
parser.add_argument('--reduced',required=False,help="Only store minimal information (LOD>4.0, no hits in metadata)")
parser.add_argument('--debug',required=False,help="Debug mode")
args = parser.parse_args()

# ASCII
X=ord('X')
Y=ord('Y')
M=ord('M')

meta = { "type": "gemma-assoc",
         "version": 1.0,
         "key-format": ">cL",
         "rec-format": "=ffff" }
log = [] # track output log
hits = [] # track hits

if args.meta:
    meta["gemma-wrapper"] = json.load(open(args.meta))

if "trait_values" in meta["gemma-wrapper"]["meta"]:
    # --- Load the traits and do some statistics
    named_values = meta["gemma-wrapper"]["meta"]["trait_values"]
    values = []
    for ind,value in named_values.items():
       values.append(value)

    a =  np.array(values)
    # print("@@@@@@",values," mean=",a.mean()," std=",a.std()," kurtosis=",stat.kurtosis(a)," skew=",stat.skew(a))
    meta["nind"] = round(a.size,4)
    meta["mean"] = round(a.mean(),4)
    meta["std"] = round(a.std(),4)
    meta["skew"] = round(stat.skew(a),4)
    meta["kurtosis"] = round(stat.kurtosis(a),4)

with lmdb.open(args.db,subdir=False,map_size=int(1e9)) as env:
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
                        elif (chr =='Y'):
                            chr = Y
                        elif (chr =='M'):
                            chr = M
                        else:
                            chr = int(chr)
                        LOD = -math.log10(float(p_lrt))
                        if not args.reduced or LOD >= SIGNIFICANT:
                            chr_c = pack('c',bytes([chr]))
                            key = pack('>cL',chr_c,int(pos))
                            val = pack('=fffff', float(af), float(beta), float(se), float(l_mle), float(p_lrt))
                            if args.debug:
                                test_chr_c,test_pos = unpack('>cL',key)
                                assert chr_c == test_chr_c
                                assert test_pos == int(pos)
                                test_chr = unpack('c',chr_c)
                                assert len(val)==20, f"Packed size is expected to be 20, but is {len(val)}"
                            res = txn.put(key, bytes(val), dupdata=False, overwrite=False)
                            if res == 0:
                                print(f"WARNING: failed to update lmdb record with key {key} -- probably a duplicate {chr}:{pos} ({test_chr_c}:{test_pos})")
                            else:
                                if not args.reduced and LOD >= SIGNIFICANT:
                                    hits.append([chr,int(pos),rs,p_lrt])

    with env.begin() as txn:
        with txn.cursor() as curs:
            # quick check and output of keys
            for key in list(txn.cursor().iternext(values=False)):
                chr,pos = unpack('>cL',key)

    meta["hits"] = hits
    meta["log"] = log
    # make it reproducible by removing variable items
    del meta["gemma-wrapper"]["time"]
    del meta["gemma-wrapper"]["user_time"]
    del meta["gemma-wrapper"]["system_time"]
    del meta["gemma-wrapper"]["wall_clock"]
    del meta["gemma-wrapper"]["ram_usage_gb"]
    del meta["gemma-wrapper"]["run_stats"]
    del meta["gemma-wrapper"]["user"]
    del meta["gemma-wrapper"]["hostname"]
    del meta["gemma-wrapper"]["gemma_command"]

    # print("HELLO: ",file=sys.stderr)
    print(meta,file=sys.stderr)
    # --- Store the metadata as a JSON record in the DB
    with env.begin(write=True) as txn:
        res = txn.put('meta'.encode(), json.dumps(meta).encode(), dupdata=False, overwrite=False)
