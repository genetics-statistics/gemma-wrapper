#!/usr/bin/env ruby
#
# Take the GWA files as input (*GWA.tar.xz) and write out to RDF. If you provide the SNP annotations their names will be plugged in.
# Otherwise it will use naming based on location (chr+pos).
#
# If you get a compatibility error in guix you may need an older Ruby. Otherwise you can do:
#
#   env GEM_PATH=tmp/ruby GEM_HOME=tmp/ruby gem install lmdb
#   env GEM_PATH=tmp/ruby ruby -e "require 'lmdb'"
#
# Pjotr Prins (c) 2025

require 'tmpdir'
require 'json'
require 'lmdb'
require 'optparse'
require 'socket'

LOD_THRESHOLD = 5.0

X='X'.ord
Y='Y'.ord
M='M'.ord

options = { show_help: false }

opts = OptionParser.new do |o|
  o.banner = "\nUsage: #{File.basename($0)} [options] filename(s)"

  o.on('-a','--anno FILEN', 'Annotation file') do |anno|
    options[:anno] = anno
    raise "Annotation input file #{anno} does not exist" if !File.exist?(anno)
  end

  o.on("--run-name NAME", "set name of run") do |b|
    options[:run_name] = b
  end

  o.on("--meta", "Output metadata only") do |b|
    options[:meta] = b
  end

  o.on("--snps", "Output snp annotation only") do |b|
    options[:snps] = b
  end

  o.on("--sort", "Sort output by significance") do |b|
    options[:sort] = b
  end

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

snps = {}

is_anno_mdb = false
snp_env = nil
db2 = nil

snpfn = options[:anno]

if snpfn
  $stderr.print "Reading #{snpfn}\n"

  if snpfn =~ /\.mdb$/
    is_anno_mdb = true
    $stderr.print("lmdb #{snpfn}...\n")
    snp_env = LMDB.new(snpfn, nosubdir: true)
    snp_db = snp_env.database("marker",create: false)
  else
    # the text file option (GEMMA annotation file)
    File.open(snpfn).each_line do |line|
      snp,pos,chr = line.split(/[\s,]+/)
      snps[chr+":"+pos] = snp
    end
  end
end

def rdf_normalize(uri)
  uri.gsub(/-/,"_")
end

USER=ENV['USER']
HOSTNAME=Socket.gethostname

used_snps = {}

get_marker_name_and_key = lambda { |chr,pos|
  location = "#{chr}:#{pos}"
  if is_anno_mdb
    chr_c =
      if chr == "X"
        X
      elsif chr == "Y"
        Y
      elsif chr == "M"
        M
      else
        chr.to_i
      end
    key = [chr_c,pos.to_i].pack("cL>")
    marker_name = snp_db[key]
    # p [chr,pos,chr_c,pos.to_i,marker_name,snp_db]
  else
    marker_name =
      if options[:anno] and snps.has_key?(location)
        snps[location]
      else
        snp = "chr#{chr}_#{pos}"
      end
    key = location
  end
  return marker_name, key
}

get_marker_info_by_key = lambda { |key|
  if is_anno_mdb
    marker_name = snp_db[key]
    chr1,pos = key.unpack('cL>')
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
    pos = pos.to_s
  else
    chr,pos = key.split(":")
    marker_name =
      if snps.has_key? key
        snps[key]
      end
  end
  return marker_name,chr,pos
}

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

      begin
        meta = JSON.parse(db['meta'])
      rescue JSON::ParserError
        next
      end

      name  = meta['gemma-wrapper']['meta']['name']
      trait = meta['gemma-wrapper']['meta']['trait']
      gwa   = meta['gemma-wrapper']['meta']['archive_GWA']
      loco  = meta['gemma-wrapper']['meta']['loco']
      xtime = meta['gemma-wrapper']['input']['time']
      run = if options[:run_name]
              options[:run_name]
            else
              run   = meta['gemma-wrapper']['input']['name']
            end
      nind  = meta['nind']
      mean  = meta['mean']
      std   = meta['std']
      skew  = meta['skew']
      kurtosis  = meta['kurtosis']
      $stderr.print("Dataset for #{name} #{trait}\n")
      result = []
      db.each do | key,value |
        chr,pos = key.unpack('cL>') # note pos is big-endian stored for easy sorting
        af,beta,se,l_mle,p_lrt = value.unpack('fffff')

        marker,location = get_marker_name_and_key.call(chr,pos)

        if not options[:snps]
          marker = "Chr#{chr}_#{pos}" if not marker
          effect = -(beta/2.0)
          minusLogP = -Math.log10(p_lrt)
          # p [p_lrt,minusLogP]
          minusLogP = 0.0 if p_lrt.nan? # not correct, but main thing is it does not show
          if minusLogP >= LOD_THRESHOLD
            rec = {chr: chr, pos: pos, snp: rdf_normalize(marker).capitalize, af: af.round(3), se: se.round(3), effect: effect.round(3), logP: minusLogP.round(2)}
            result.push rec
          end
        end
        used_snps[location] = true
      end
      env.close
      if not options[:snps] # output all triples
        if options[:sort]
          $stderr.print("Sorting...\n")
          result = result.sort_by { |rec| rec[:logP] }.reverse
        end

        @prefix = "GEMMAMapped"
        hash = gwa[32..39]
        postfix = rdf_normalize(gwa[41..-8])+"_"+hash
        s_loco = (loco ? "LOCO" : "")
        id = "gn:#{@prefix}_#{run}_#{s_loco}_#{postfix}"
        print """#{id} a gnt:mappedTrait;
      rdfs:label \"GEMMA #{name} trait #{trait} mapped with LOCO (defaults)\";
      gnt:trait gn:publishXRef_#{trait};
      gnt:loco #{loco};
      gnt:run gn:#{run};
      gnt:time \"#{xtime}\";
      gnt:belongsToGroup gn:setBxd;
      gnt:name \"#{name}\";
      gnt:traitId \"#{trait}\";
      gnt:nind #{nind};
      gnt:mean #{mean};
      gnt:std #{std};
      gnt:skew #{skew};
      gnt:kurtosis #{kurtosis};
      skos:altLabel \"BXD_#{trait}\";
      gnt:filename \"#{File.basename(fn)}\";
      gnt:hostname \"#{HOSTNAME}\";
      gnt:user \"#{USER}\".
"""

        first = true
        result.each do |rec|
          # we always show the highest hit
          locus = rec[:snp]
          locus="unknown" if locus == "?"
          # for the rest of the hits make sure they are significant and have a snp id:
          if not first
            break if rec[:logP] < LOD_THRESHOLD # superfluous
            next if locus == "unknown"
          end
          # rdfs:label \"Mapped locus #{locus} for #{name} #{trait}\";
          # FIXME locus.capitalize
          print """gn:#{locus}_#{postfix} a gnt:mappedLocus;
      gnt:mappedSnp #{id};
      gnt:locus gn:#{locus.capitalize};
      gnt:lodScore #{rec[:logP].round(1)};
      gnt:af #{rec[:af]};
      gnt:effect #{rec[:effect]}.
"""
          first = false
        end
      end
    end
  end
end # tmpdir

# p used_snps

if options[:snps] # output SNP annotation only
  used_snps.each_key do | key |
    name,chr,pos = get_marker_info_by_key.call(key)
    if name
      print """gn:#{rdf_normalize(name)} a gnt:marker;
                 rdfs:label \"#{name}\";
                 gnt:chr  \"#{chr}\";
                 gnt:pos  #{pos}.
"""
    end
  end
end
