#! /bin/env sh
#
# Take a list of IDs and fetch SNPs


export TMPDIR=./tmp
export RDF=pan-qtl.rdf
./bin/gemma2rdf.rb --header > $RDF

for id in `head -1 pan-ids.txt` ; do
    echo Precomputing $0 for $id
    curl -G http://sparql-test.genenetwork.org/sparql -H "Accept: application/json; charset=utf-8" --data-urlencode query="
PREFIX dct: <http://purl.org/dc/terms/>
PREFIX gn: <http://genenetwork.org/id/>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX gnc: <http://genenetwork.org/category/>
PREFIX gnt: <http://genenetwork.org/term/>
PREFIX sdmx-measure: <http://purl.org/linked-data/sdmx/2009/measure#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX qb: <http://purl.org/linked-data/cube#>
PREFIX xkos: <http://rdf-vocabulary.ddialliance.org/xkos#>
PREFIX pubmed: <http://rdf.ncbi.nlm.nih.gov/pubmed/>

SELECT * WHERE { ?traitid a gnt:mappedTrait ; gnt:traitId ?trait ; gnt:kurtosis ?k . } LIMIT 3
"

done
