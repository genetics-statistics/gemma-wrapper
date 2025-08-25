#! /bin/env sh

export TMPDIR=./tmp
curl http://127.0.0.1:8092/dataset/bxd-publish/list > bxd-publish.json
jq ".[] | .Id" < bxd-publish.json > ids.txt

for id in `cat ids.txt` ; do
  echo Precomputing $id
  if [ ! -e tmp/trait-BXDPublish-$id-gemma-GWA-hk.tar.xz ] ; then
      curl http://127.0.0.1:8092/dataset/bxd-publish/values/$id.json > pheno.json
      ./bin/gn-pheno-to-gemma.rb --phenotypes pheno.json --geno-json BXD.geno.json > BXD_pheno.txt
      # ./bin/gemma-wrapper --json --lmdb --geno-json BXD.geno.json --phenotypes pheno.json --population BXD --name BXDPublish --trait $id --loco --input K.json -- -g BXD.geno.txt -a BXD.8_snps.txt -lmm 9 -maf 0.1 -n 2 -debug > GWA.json
      gemma -g BXD.geno.txt -p BXD_pheno.txt -a BXD.8_snps.txt -n 2 -lm 2 -o tmp/trait-BXDPublish-$id-gemma-GWA-hk
  fi
done
