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
does a pass-through of all standard GEMMA invocation switches.

## Installation

Prerequisites are

* A recent version of [GEMMA](https://github.com/genetics-statistics/GEMMA)
* Standard [Ruby >2.0 ](https://www.ruby-lang.org/en/) which comes on
  almost all Linux systems

Fetch a [release](https://github.com/genetics-statistics/gemma-K-handler/releases) of
[gemma-k-handler](https://github.com/genetics-statistics/gemma-K-handler)

Unpack it and run the tool as

    ./bin/gemma-K-handler --help

## Usage

## Copyright

Copyright (c) 2017 Pjotr Prins. See LICENSE.txt for further details.
