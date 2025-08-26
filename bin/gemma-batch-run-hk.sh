#! /bin/env sh

export TMPDIR=./tmp
export GN_GUILE=http://127.0.0.1:8091
export RDF=gemma-GWA-hk.ttl
curl $GN_GUILE/dataset/bxd-publish/list > bxd-publish.json
jq ".[] | .Id" < bxd-publish.json |head -10 > ids.txt
./bin/gemma2rdf.rb --header > $RDF

for id in `cat ids.txt` ; do
    echo Precomputing $0 for $id
    traitfn=trait-BXDPublish-$id-gemma-GWA-hk.assoc.txt
    if [ ! -e $TMPDIR/$trainfn ] ; then
        curl $GN_GUILE/dataset/bxd-publish/values/$id.json > pheno.json
        ./bin/gn-pheno-to-gemma.rb --phenotypes pheno.json --geno-json BXD.geno.json > BXD_pheno.txt
        gemma -g BXD.geno.txt -p BXD_pheno.txt -a BXD.8_snps.txt -n 2 -lm 2 -outdir $TMPDIR -o $traitfn
    fi
    ./bin/gemma2rdf.rb $TMPDIR/$traitfn >> $RDF

done
