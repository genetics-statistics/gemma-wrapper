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

qtls = QTL::QRanges.new("10002","test")
ARGV.each do |fn|
  CSV.foreach(fn,headers: true, col_sep: "\t") do |hit|
    lod = hit["lod"].to_f
    if lod > 5.0
      qlocus = QTL::QLocus.new(hit["nodeid"],hit["chr"],hit["pos"].to_f/10**6,hit["af"].to_f,lod)
      qtls.add_locus(qlocus)
    end
  end
end
qtls.rebin
qtls.pangenome_filter

print qtls
