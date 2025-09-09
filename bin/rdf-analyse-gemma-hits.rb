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

require 'gnrdf'
require 'qtlrange'

options = { show_help: false }

opts = OptionParser.new do |o|
  o.banner = "\nUsage: #{File.basename($0)} [options] assoc filename(s)"

  o.on_tail('--header', 'Write header') do
    options[:header] = true
  end

  o.on('-o','--output TYPE', 'Output TEXT (default) or RDF') do |type|
    options[:output] = type.to_sym
    options[:rdf] = true if type == "RDF"
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

=begin
Store trait info + snps, e.g 'p loco':
 "10008"=>
  {:loco=>true,
   :traitId=>"10008",
   :snps=>
    [#<RDF::URI:0xd73c URI:http://genenetwork.org/id/Rsm10000011798_BXDPublish_10008_gemma_GWA_7c00f36d>,
     #<RDF::URI:0xd7a0 URI:http://genenetwork.org/id/Rs48427909_BXDPublish_10008_gemma_GWA_7c00f36d>,
     #<RDF::URI:0xd804 URI:http://genenetwork.org/id/Rs31915734_BXDPublish_10008_gemma_GWA_7c00f36d>],
   :id=>#<RDF::URI:0xd6b0 URI:http://genenetwork.org/id/GEMMAMapped_LOCO_BXDPublish_10008_gemma_GWA_7c00f36d>}}
=end

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
    $stderr.print [traitid,combined.size,difference.size],"\n"
    # let's try to define ranges
    # we need to revert to using SNP unique IDs now

    gemma_snps = loco[traitid][:snps]
    hk_snps = hk[traitid][:snps]

    if difference.size > 0
      results = {}
      [[combined,"combined"],[hk_snps,"HK"],[gemma_snps,"LOCO"]].each do |cmd|
        set,setname = cmd
        qtls = QTL::QRanges.new(traitid,setname)
        set.each do | snp |
          base_snp_id = snp
          snp_info = snps[base_snp_id.to_s]
          if snp_info == nil
            # redirect to the actual snp id
            base_snp_id = locus[snp].to_s
            snp_info = snps[base_snp_id.to_s]
          end
          snp_uri = snp_info["snp"]
          chr = snp_info["chr"]
          pos = snp_info["mb"].to_f
          snp_lod = lod[snp]
          snp_lod = snp_lod.to_f if snp_lod != nil
          qlocus = QTL::QLocus.new(snp.to_s,chr,pos,snp_lod)
          qtls.add_locus(qlocus)
        end
        results[setname] = qtls
      end
      id = gnid(loco[traitid][:id])
      results["LOCO"].print_rdf(id) if options[:rdf]
      qtl_diff(id,results["HK"],results["LOCO"],options[:rdf])
    end
  else
    $stderr.print "WARNING: no HK counterpart for #{traitid}\n"
  end
end
