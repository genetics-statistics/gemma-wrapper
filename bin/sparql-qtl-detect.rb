#! /bin/env ruby
#
# Detects the QTLs for a trait from a list of markers/SNPs with their LOD scores and outputs the annotated QTL as RDF

require 'csv'
require 'optparse'
require 'pp'

basepath = File.dirname(File.dirname(__FILE__))
$: << File.join(basepath,'lib')

require 'gnrdf'
require 'qtlrange'

options = { show_help: false }

opts = OptionParser.new do |o|
  o.banner = "\nUsage: #{File.basename($0)} [options] filename(s)"

  o.on('-o','--output TYPE', 'Output TEXT (default) or RDF') do |type|
    options[:output] = type.to_sym
  end

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

OUTPUT_RDF = options[:output] == :RDF

if options[:header]
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

ARGV.each do |fn|
  trait_id = nil
  trait = fn.split(".")[0]
  qtls = QTL::QRanges.new(trait,"test")
  CSV.foreach(fn,headers: true, col_sep: "\t") do |hit|
    trait_id = hit["traitid"] if not trait_id
    lod = hit["lod"].to_f
    if lod > 5.0 # set for pangenome input
      qlocus = QTL::QLocus.new(hit["snp"],hit["chr"],hit["pos"].to_f/10**6,hit["af"].to_f,lod)
      qtls.add_locus(qlocus)
    end
  end
  qtls.rebin
  qtls.pangenome_filter
  if OUTPUT_RDF
    qtls.print_rdf trait,trait_id
  else
    print qtls
  end
end
