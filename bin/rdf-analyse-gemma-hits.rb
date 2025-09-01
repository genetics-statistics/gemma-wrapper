#!/usr/bin/env ruby
#
# Parse SNP RDF and compare between two datasets
#
# Pjotr Prins (c) 2025

require 'tmpdir'
require 'rdf'
require 'rdf/turtle'
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

if options[:show_help] or ARGV.size == 0
  print opts
  # print USAGE
  exit 1
end

GN = RDF::Vocabulary.new("http://genenetwork.org/id/")
GNT = RDF::Vocabulary.new("http://genenetwork.org/term/")

ARGV.each do | fn |
  $stderr.print "Parsing #{fn}...\n"
  require 'rdf/ntriples'
  graph = RDF::Graph.load(fn)

  datasets = graph.query(RDF::Query.new {
                           pattern [:dataset, RDF.type, GNT.mappedTrait]
                         })

  datasets.each { |trait|
    p "-------"
    p trait.dataset
    snps = graph.query(RDF::Query.new {
                         pattern [ :snp, GNT.mappedSnp, trait.dataset ]
                       })
    p snps
  }

end
