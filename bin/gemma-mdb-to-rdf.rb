#!/usr/bin/env ruby
#
# To run this you'll need to do:
#
#   env GEM_PATH=tmp/ruby GEM_HOME=tmp/ruby gem install lmdb
#   env GEM_PATH=tmp/ruby GEM_HOME=tmp/ruby ruby -e "require 'lmdb'"
#
# until I fixed the package.
#
# Pjotr Prins (c) 2025

require 'tmpdir'
require 'json'
require 'lmdb'
require 'optparse'

options = { show_help: false }

opts = OptionParser.new do |o|
  o.banner = "\nUsage: #{File.basename($0)} [options] filename(s)"

  o.on('-a','--anno FILEN', 'Annotation file') do |anno|
    options[:anno] = anno
    raise "Annotation input file #{anno} does not exist" if !File.exist?(anno)
  end

  o.on("--meta", "Output metadata only") do |b|
    options[:meta] = b
  end

  o.on("--sort", "Sort output by significance") do |b|
    options[:sort] = b
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

snps = {}

if options[:anno]
  File.open(options[:anno]).each_line do |line|
    snp,pos,chr = line.split(/[\s,]+/)
    snps[chr+":"+pos] = snp
  end
end

def rdf_normalize(uri)
  uri.gsub(/-/,"_")
end

ARGV.each do |fn|
  Dir.mktmpdir do |tmpdir|
    $stderr.print("Parsing #{fn}...\n")
    if fn =~ /xz$/
      $stderr.print `tar xvJf #{fn} -C #{tmpdir} > /dev/null`
    else
      raise "Expected xz tarball!"
    end
    Dir.glob(tmpdir+"/*.mdb").each do |mdb|
      $stderr.print("lmdb #{mdb}...\n")
      env = LMDB.new(mdb, nosubdir: true)
      maindb = env.database
      db = env.database(File.basename(mdb))

      if options[:meta]
        print db['meta']
        env.close
        exit 0
      end

      meta = JSON.parse(db['meta'])

      name = meta['gemma-wrapper']['meta']['name']
      trait = meta['gemma-wrapper']['meta']['trait']
      gwa = meta['gemma-wrapper']['meta']['archive_GWA']
      loco = meta['gemma-wrapper']['meta']['loco']
      xtime = meta['gemma-wrapper']['input']['time']
      $stderr.print("Dataset for #{name} #{trait}\n")
      result = []
      db.each do | key,value |
        chr,pos = key.unpack('cL>') # note pos is big-endian stored for easy sorting
        af,beta,se,l_mle,p_lrt = value.unpack('fffff')
        snp = "?"
        location = "#{chr}:#{pos}"
        if options[:anno] and snps.has_key?(location)
          snp = snps[location]
        end
        effect = -(beta/2.0)
        minusLogP = -Math.log10(p_lrt)
        # p [p_lrt,minusLogP]
        minusLogP = 0.0 if p_lrt.nan?
        rec = {chr: chr, pos: pos, snp: rdf_normalize(snp).capitalize, af: af.round(3), se: se.round(3), effect: effect.round(3), logP: minusLogP.round(2)}
        result.push rec
      end
      env.close
      if options[:sort]
        $stderr.print("Sorting...\n")
        result = result.sort_by { |rec| rec[:logP] }.reverse
      end

      print """
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
"""

      prefix = "GEMMAMapped"
      hash = gwa[32..39]
      postfix = rdf_normalize(gwa[41..-8])+"_"+hash
      s_loco = (loco ? "LOCO" : "")
      id = "gn:#{prefix}_#{s_loco}_#{postfix}"
      print """#{id} a gnt:mappedTrait;
      rdfs:label \"GEMMA #{name} trait #{trait} mapped with LOCO (defaults)\";
      gnt:trait gn:publishXRef_#{trait};
      gnt:loco #{loco};
      gnt:time \"#{xtime}\";
      gnt:belongsToGroup gn:setBxd;
      gnt:name \"#{name}\";
      gnt:traitId \"#{trait}\";
      skos:altLabel \"BXD_#{trait}\".
"""

      first = true
      result.each do |rec|
        # we always show the highest hit
        locus = rec[:snp]
        locus="unknown" if locus == "?"
        # for the rest of the hits make sure they are significant and have a snp id:
        if not first
          break if rec[:logP] < 4.0
          next if locus == "unknown"
        end
      # rdfs:label \"Mapped locus #{locus} for #{name} #{trait}\";
        print """gn:#{locus}_#{postfix} a gnt:mappedLocus;
      gnt:mappedSnp #{id};
      gnt:locus gn:#{locus};
      gnt:lodScore #{rec[:logP]};
      gnt:af #{rec[:af]};
      gnt:effect #{rec[:effect]}.
"""
        first = false
      end
    end
  end
end # tmpdir
