#!/usr/bin/env ruby
#
# Parse SNP RDF and compare between two datasets
#
# Pjotr Prins (c) 2025

require 'csv'
require 'tmpdir'
require 'rdf'
require 'rdf/turtle'
require 'rdf/raptor'
require 'optparse'
require 'pp'

basepath = File.dirname(File.dirname(__FILE__))
$: << File.join(basepath,'lib')

require 'qtlrange'

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

snps = {}
CSV.foreach("snps.txt",headers: true, col_sep: "\t") do |row|
  snps[row["snp"]] = row
end

GN = RDF::Vocabulary.new("http://genenetwork.org/id/")
GNT = RDF::Vocabulary.new("http://genenetwork.org/term/")

# graph = RDF::Graph.new
traits = {}
locus = {}
lod = {}
ARGV.each do | fn |
  $stderr.print "Parsing #{fn}...\n"
  reader = RDF::Reader.open(fn)
  reader.each_statement do |statement|
    # p statement.inspect
    # these are stored by run-trait:
    subject = statement.subject
    traits[subject] = {} if statement.object == GNT.mappedTrait
    traits[subject][:traitId] = statement.object.to_s if statement.predicate == GNT.traitId
    traits[subject][:loco] = true if statement.predicate == GNT.loco
    traits[subject][:hk] = true if statement.predicate == GNT.gemmaHk
    # note we assume SNPs come after! Store by run-trait-SNP id:
    if statement.predicate == GNT.mappedSnp
      traitid = statement.object
      traits[traitid][:snps] ||= []
      traits[traitid][:snps].push statement.subject
    end
    locus[statement.subject] = statement.object if statement.predicate == GNT.locus
    lod[statement.subject] = statement.object if statement.predicate == GNT.lodScore # locus and lod share identifier
  end
end

$stderr.print "# traits is #{traits.size}\n"

loco = {}
hk = {}
traits.each do |k,v|
  traitid = v[:traitId]
  v[:id] = k
  loco[traitid] = v if v[:loco]
  hk[traitid] = v if v[:hk]
end

$stderr.print "# loco LMM is #{loco.size}\n"
$stderr.print "# HK is #{hk.size}\n"

# Now we have two sets of traits and we will walk every GEMMA set, and see if it matches HK

loco.each do | traitid, rec |
  # Walk all traidid and see if we have an HK counterpart
  if hk[traitid]
    # We create two sets out of the SNPs, make sure to transform to locus names for comparison
    gemma_snps = loco[traitid][:snps].map { |snp| locus[snp] }
    hk_snps = hk[traitid][:snps].map { |snp| locus[snp] }
    gemma_set = Set.new(gemma_snps)
    hk_set = Set.new(hk_snps)
    combined = gemma_set + hk_set
    difference = gemma_set - hk_set
    p [traitid,combined.size,difference.size]
    # let's try to define ranges
    # we need lod scores

    #gemma_snps_lod = loco[traitid][:snps].map { |snp| lod[snp].to_f }
    #hk_snps_lod = loco[traitid][:snps].map { |snp| lod[snp].to_f }

    if difference.size > 0
      [[combined,"combined"],[hk_set,"HK"],[gemma_set,"LOCO"]].each do |cmd|
        set,setname = cmd
        qtls = QTL::QRanges.new(traitid,setname)
        set.each do | snp |
          snp_info = snps[snp.to_s]
          snp_uri = snp_info["snp"]
          chr = snp_info["chr"]
          pos = snp_info["mb"].to_f
          qlocus = QTL::QLocus.new(snp_uri,chr,pos)
          qtls.add_locus(qlocus)
        end
        p qtls
      end
    end
  else
    $stderr.print "WARNING: no HK counterpart for #{traitid}\n"
  end
end
