#!/usr/bin/env ruby
#
# This is a simplified and flexible version of ./gemma-mdb-rdf.rb -- without lmdb
#
#     ./bin/gemma2rdf.rb output/trait-BXDPublish-*-gemma-GWA-hk.assoc.txt
#
# will write RDF to stdout
#
# Pjotr Prins (c) 2025

require 'tmpdir'
require 'json'
require 'lmdb'
require 'optparse'

options = { show_help: false }

opts = OptionParser.new do |o|
  o.banner = "\nUsage: #{File.basename($0)} [options] assoc filename(s)"

  o.on_tail('--header', 'Write header') do
    options[:header] = true
  end

  o.separator ""

  o.on_tail('-h', '--help', 'display this help and exit') do
    options[:show_help] = true
  end
end

opts.parse!(ARGV)

if options[:show_help]
  print opts
  # print USAGE
  exit 1
end

if options[:header]

# Other prefixes used in our store:
# @prefix pubmed: <http://rdf.ncbi.nlm.nih.gov/pubmed/> .
# @prefix qb: <http://purl.org/linked-data/cube#> .
# @prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
# @prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
# @prefix sdmx-measure: <http://purl.org/linked-data/sdmx/2009/measure#> .
# @prefix skos: <http://www.w3.org/2004/02/skos/core#> .
# @prefix xkos: <http://rdf-vocabulary.ddialliance.org/xkos#> .
# @prefix xsd: <http://www.w3.org/2001/XMLSchema#> .

      print """
@prefix dct: <http://purl.org/dc/terms/> .
@prefix gn: <http://genenetwork.org/id/> .
@prefix gnc: <http://genenetwork.org/category/> .
@prefix gnt: <http://genenetwork.org/term/> .
@prefix owl: <http://www.w3.org/2002/07/owl#> .
@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
@prefix skos: <http://www.w3.org/2004/02/skos/core#> .
"""
end

def rdf_normalize(uri)
  uri.gsub(/\W/,"_")
end

TIMEX = Time.now.to_s[..30]

ARGV.each do |fn|
  $stderr.print("Parsing #{fn}\n")
  header = nil
  recs = []
  File.open(fn).each_line do |line|
    # first read the records because we only want to output the top hits
    rec = line.strip.split(/\t/)
    if rec[0] == "chr"
      header = rec
    else
      minusLogP = -Math.log10(rec[-1].to_f)
      rec.push minusLogP
      recs.push rec
    end
  end
  header.push "LOD"
  # p [header,recs.size]
  col = header.size - 1
  sorted = recs.sort_by { |rec| rec[-1].to_f }.reverse
  name = "BXDPublish"
  fn =~ /(\d+)/
  trait = $1
  id = rdf_normalize("HK_#{File.basename(fn)}")
  print """
gn:#{id} a gnt:mappedTrait;
        rdfs:label \"GEMMA_BXDPublish #{fn} trait HK mapped\";
        gnt:GEMMA_HK true;
        gnt:belongsToGroup gn:setBxd;
        gnt:trait gn:publishXRef_#{trait};
        gnt:time \"#{TIMEX}\";
        gnt:name \"#{name}\";
        gnt:traitId \"#{trait}\";
        skos:altLabel \"BXD_#{trait}\".
"""
  sorted.each do | rec |
    lod = rec[-1]
    if lod > 3.0
      chr,snpname = rec
      snp = "#{snpname}_#{id}"
      # FIXME snpname.capitalize
      print """
gn:#{snp} a gnt:mappedLocus;
       gnt:mappedSnp gn:#{snp} ;
       gnt:locus gn:#{snpname.capitalize} ;
       gnt:lodScore #{lod.round(2)} .
"""
    end
  end
end
