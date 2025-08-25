#! /bin/env sh

export TMPDIR=./tmp
curl http://127.0.0.1:8092/dataset/bxd-publish/list > bxd-publish.json
jq ".[] | .Id" < bxd-publish.json > ids.txt
./bin/gemma-wrapper --force --json --loco -- -g BXD.geno.txt -p BXD_pheno.txt -a BXD.8_snps.txt -n 2 -gk > K.json

for id in `cat ids.txt` ; do
  echo Precomputing $id 
  if [ ! -e tmp/*-BXDPublish-$id-gemma-GWA.tar.xz ] ; then
    curl http://127.0.0.1:8092/dataset/bxd-publish/values/$id.json > pheno.json
    ./bin/gemma-wrapper --json --lmdb --geno-json BXD.geno.json --phenotypes pheno.json --population BXD --name BXDPublish --trait $id --loco --input K.json -- -g BXD.geno.txt -a BXD.8_snps.txt -lmm 9 -maf 0.1 -n 2 -debug > GWA.json
  fi
done
