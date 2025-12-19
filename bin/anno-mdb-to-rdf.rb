#!/usr/bin/env ruby
#
# Read anno mdb (created by anno2mdb.rb) and write to RDF.
#
# Pjotr Prins (c) 2025

require 'tmpdir'
require 'json'
require 'lmdb'
require 'optparse'
require 'socket'

LOD_THRESHOLD = 5.0
CHRPOS_PACK="S>L>L>" # chr, pos, line. L is uint32, S is uint16 - total 64bit

X='X'.ord
Y='Y'.ord
M='M'.ord

options = { show_help: false }

opts = OptionParser.new do |o|
  o.banner = "\nUsage: #{File.basename($0)} [options] filename"

  o.on_tail('--header', 'Write header') do
    options[:header] = true
  end

  # o.on("-v", "--verbose", "Run verbosely") do |v|
  #   options[:verbose] = true
  # end

  # o.on("-d", "--debug", "Show debug messages and keep intermediate output") do |v|
  #   options[:debug] = true
  # end

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
  uri.gsub(/[-<+>]/,"_")
end

annofn = ARGV.shift
anno_marker_tab = nil

$stderr.print "Reading #{annofn}\n"
marker_env = LMDB.new(annofn, nosubdir: true)
begin
  anno_marker_tab = marker_env.database("marker", create: false)
  marker_info = marker_env.database("info",create: false)
  raise "Metadata missing in anno mdb file" if not marker_info["meta"]
  meta = JSON.parse(marker_info["meta"])
  raise "Not an anno file #{meta}" if meta["type"] != "gemma-anno"
  raise "Incompatible key format" if meta["key-format"] != CHRPOS_PACK
rescue
  raise "Problem reading annotation file #{annofn}!"
end

anno_marker_tab.each do |key,marker|
  # chr1,pos = key.unpack('cL>')
  chr1,pos,num = key.unpack(CHRPOS_PACK)
  # p [key,marker,chr1,pos,num]
  chr =
    if chr1 == X
      "X"
    elsif chr1 == Y
      "Y"
    elsif chr1 == M
      "M"
    else
      chr1
    end
  # pos = pos.to_f/1_000_000 - let's not adjust
  id = rdf_normalize(marker)
  print """gn:#{id} a gnt:locus;
      rdfs:label \"#{marker}\";
      gnt:chr \"#{chr}\" ;
      gnt:pos #{pos} .
"""
end
