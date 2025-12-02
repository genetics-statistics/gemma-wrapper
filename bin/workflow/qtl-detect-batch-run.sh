#! /bin/env sh
#
# Take a list of IDs and fetch SNPs
#


export TMPDIR=./tmp
export RDF=pan-qtl.rdf

for id in `cat pan-ids-sorted.txt` ; do
    echo Precomputing $0 for $id
    if [ -e $id.hits.txt ] ; then
        continue
    fi
    echo "gnt:traitId $id ."

    curl -G http://sparql-test.genenetwork.org/sparql -H "Accept: text/tab-separated-values; charset=utf-8" --data-urlencode query="
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

SELECT ?traitid ?lod ?af ?snp ?chr ?pos FROM <http://pan-test.genenetwork.org> WHERE {
?traitid a gnt:mappedTrait;
         gnt:run gn:test ;
         gnt:traitId \"$id\" ;
         gnt:traitId ?trait .
?snp gnt:mappedSnp ?traitid ;
        gnt:locus ?locus ;
        gnt:lodScore ?lod ;
        gnt:af ?af .
FILTER(?lod >= 5.0) .
?locus rdfs:label ?marker ;
         gnt:chr ?chr ;
         gnt:pos ?pos .
FILTER (contains(?marker,\"Marker\") && ?pos > 1000) # FIXME: this is to avoid duplicates
} ORDER BY DESC(?lod)
" > $id.hits.txt

  ../../bin/sparql-qtl-detect.rb --header $id.hits.txt -o RDF >> $RDF

done
