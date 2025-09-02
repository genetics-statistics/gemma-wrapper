#!/usr/bin/env ruby
#
# Parse SNP RDF and compare between two datasets
#
# Pjotr Prins (c) 2025

require 'tmpdir'
require 'rdf'
require 'rdf/turtle'
require 'rdf/raptor'

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

# graph = RDF::Graph.new
traits = {}
snps = {}
ARGV.each do | fn |
  $stderr.print "Parsing #{fn}...\n"
  reader = RDF::Reader.open(fn)
  reader.each_statement do |statement|
    # p statement.inspect
    subject = statement.subject
    traits[subject] = {} if statement.object == GNT.mappedTrait
    traits[subject][:traitId] = statement.object.to_s if statement.predicate == GNT.traitId
    traits[subject][:loco] = statement.object.to_s if statement.predicate == GNT.loco
    # note we assume SNPs come after!
    if statement.predicate == GNT.mappedSnp
      traitid = statement.object
      traits[traitid][:snps] ||= []
      traits[traitid][:snps].push statement.subject
    end
    snps[statement.subject] = statement.object if statement.predicate == GNT.locus
  end
end

# p traits
# p traits.size
# p snps
