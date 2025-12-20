#! /bin/env ruby
#
# Converts the 'epoch' spreadsheet to RDF

require 'csv'
require 'optparse'
require 'pp'

basepath = File.dirname(File.dirname(__FILE__))
$: << File.join(basepath,'lib')

require 'gnrdf'

options = { show_help: false }

opts = OptionParser.new do |o|
  o.banner = "\nUsage: #{File.basename($0)} [options] filename(s)"

  # o.on('-o','--output TYPE', 'Output TEXT (default) or RDF') do |type|
  #   options[:output] = type.to_sym
  # end

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

def rdf_normalize(uri)
  uri.gsub(/[-<+>\/]/,"_")
end

def rdf_expand_str(meta,key,predicate=nil,value=nil)
  if meta[key] and meta[key] != ""
    predicate = rdf_normalize(key.to_s) if not predicate
    object = meta[key]
    object = value.call(object) if value
    if predicate =~ /:/
      print "                                #{predicate} \"#{object}\" ;\n"
    else
      print "                                gnt:#{predicate} \"#{object}\" ;\n"
    end
  end
end

def rdf_expand(meta,key,predicate=nil,value=nil)
  if meta[key] and meta[key] != ""
    predicate = rdf_normalize(key.to_s) if not predicate
    object = meta[key]
    object = value.call(object) if value
    print "                                gnt:#{predicate} #{object} ;\n"
  end
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

alt = {}
info = {}
ARGV.each do |fn|
  CSV.foreach(fn,headers: true, col_sep: "\t") do |tsv|
    gn_info = tsv["Strain GeneNetwork name"]
    if gn_info
      gn_info.strip!
      w1,w2 = gn_info.split
      id = w1.strip
      if w2
        next if w1 == "Renamed"
        w2 =~ /(\w+)/
        alt[$1] = w1
      end
      meta = {
        expanded_name: tsv["Strain expanded name"],
        start_year: tsv["Year breeding started"].to_i,
        jax_start: tsv["Date production began at JAX"],
        jax_stock: tsv["JAX Stock No (RRID)"],
        uthsc_live_2023: tsv["Live at UTHSC January 2023"],
        availability_2021: tsv["Availability (as of October 2021)"],
        availability_2023: tsv["Availability (as of August 2023)"],
        epoch: tsv["Epoch"],
        method: tsv["Method of derivation"],
        birth_seq_ind: tsv["Date of birth of sequenced individual"],
        age_seq_ind: tsv["Age of individual sequenced (days)"],
        gen_seq: tsv["Generation at sequencing"],
        seq_depth: tsv["Sequencing depth of individual sequenced"],
        is_illumina_13377: tsv["Strain genotyped on Illumina 13377 SNP array"],
        is_affy_600k: tsv["Strain genotyped on Affymetrix MouseDiversityArray600K"],
        gen_aff_600k: tsv["Strain generation at Affymetrix MouseDiversityArray600K"],
        gen_muga: tsv["Strain generation at MUGA array genotyping"],
        gen_megamuga: tsv["Strain generation at MegaMUGA array genotyping"],
        gen_gigamuga: tsv["Strain generation at GigaMUGA array genotyping"],
        extinct: tsv["Generation went extinct (if known)"],
        notes: tsv["Backcross and RIX-derived breeding notes"],
        has_genotypes: tsv["Any genotypes available"],
        m_origin: tsv["Mitocondrial origin (if known)"],
        y_origin: tsv["Y-chromosome origin (if known)"],
      }
      info[id] = meta
      # p tsv
      # p [id,info[id]]
      print
          sample = rdf_normalize(id.capitalize)
          print """gn:#{sample}
                                dct:description \"#{meta[:expanded_name]}\" ;
                                gnt:epoch #{meta[:epoch]} ;
                                gnt:availability \"#{meta[:availability_2023]}\" ;
"""
          rdf_expand_str(meta,:method)
          rdf_expand_str(meta,:m_origin,"M_origin")
          rdf_expand_str(meta,:y_origin,"Y_origin")
          rdf_expand_str(meta,:jax_stock,predicate="JAX")
          rdf_expand(meta,:start_year)
          rdf_expand(meta,:age_seq_ind)
          rdf_expand_str(meta,:birth_seq_ind)
          rdf_expand_str(meta,:availability_2023)
          rdf_expand_str(meta,:extinct)
          rdf_expand_str(meta,:notes,"rdfs:comment")
          rdf_expand(meta,:has_genotypes,nil,lambda { |v| v ? "true" : "false" })
          print """                                rdfs:label \"#{id}\" .
"""
    end
  end
end
alt.each do | k,v |
  print """
  gn:#{rdf_normalize(k.capitalize)} owl:sameAs gn:#{rdf_normalize(v.capitalize)} ."""
end
