# GEMMA wrapper caches K between runs with LOCO support

![Genetic associations identified in CFW mice using GEMMA (Parker et al,
Nat. Genet., 2016)](cfw.gif)

## Introduction

GEMMA is a software toolkit for fast application of linear mixed
models (LMMs) and related models to genome-wide association studies
(GWAS) and other large-scale data sets.

This repository contains gemma-k-handler, essentially a wrapper of
GEMMA that provides support for caching the kinship or relatedness
matrix (K) with the option of full leave-one-chromosome-out genome
scans (LOCO).

gemma-k-handler requires a recent version of GEMMA and essentially
does a pass-through of all standard GEMMA invocation switches. On
return gemma-k-handler can return a JSON object (--json) which is
useful for web-services.

Note that this a work in progress (WIP). What is described below
should work.

## Installation

Prerequisites are

* A recent version of [GEMMA](https://github.com/genetics-statistics/GEMMA)
* Standard [Ruby >2.0 ](https://www.ruby-lang.org/en/) which comes on
  almost all Linux systems

Fetch a [release](https://github.com/genetics-statistics/gemma-K-handler/releases) of
[gemma-k-handler](https://github.com/genetics-statistics/gemma-K-handler)

Unpack it and run the tool as

    ./bin/gemma-k-handler --help

## Usage

gemma-k-handler picks up GEMMA from the PATH. To override that behaviour
use the GEMMA_COMMAND environment variable, e.g.

    env GEMMA_COMMAND=~/opt/gemma/bin/gemma ./bin/gemma-K-handler --help

to pass switches to GEMMA put them after '--' e.g.

    gemma-k-handler -v -- -h

prints the GEMMA help

## Caching of K

To compute K

    gemma-k-handler -- \
    -g test/data/input/BXD_geno.txt.gz \
    -p test/data/input/BXD_pheno.txt \
    -gk \
    -debug

Run it twice to see

    /tmp/3079151e14b219c3b243b673d88001c1675168b4.log.txt gemma-k-handler CACHE HIT!

gemma-k-handler computes the unique HASH value over the command
line switches passed into GEMMA as well as the contents of the files
passed in (here the genotype and phenotype files).

You can also get JSON output on STDOUT by providing the --json switch

    gemma-k-handler --json -- \
    -g test/data/input/BXD_geno.txt.gz \
    -p test/data/input/BXD_pheno.txt \
    -gk \
    -debug

prints out something that can be parsed with a calling program

```json
{"warnings":[],"errno":0,"gemma_command":"../gemma/bin/gemma -g test/data/input/BXD_geno.txt.gz -p test/data/input/BXD_pheno.txt -gk -debug -o 18ce786ab92064a7ee38a7422e7838abf91f5eb0 -outdir /tmp","type":"K","log":"/tmp/18ce786ab92064a7ee38a7422e7838abf91f5eb0.log.txt","K":"/tmp/18ce786ab92064a7ee38a7422e7838abf91f5eb0.cXX.txt"}
```

Note that GEMMA's -o (output) and --outdir switches should not be
used. gemma-k-handler stores the cached matrices in TMPDIR by
default. If you want something else provide a --cache-dir, e.g.

    gemma-k-handler --cache-dir ~/.gemma-cache -- \
    -g test/data/input/BXD_geno.txt.gz \
    -p test/data/input/BXD_pheno.txt \
    -gk \
    -debug

will write the new matrix in ~/.gemma-cache.

### LOCO

(not yet implemented)

Recent versions of GEMMA have LOCO support for a single chromosome using
the -loco switch. To loop all chromosomes gemma-k-handler do

    gemma-k-handler --loco-all -- \
        -g ../example/BXD_geno.txt.gz \
        -gk \
        -debug

## Copyright

Copyright (c) 2017 Pjotr Prins. See LICENSE.txt for further details.
