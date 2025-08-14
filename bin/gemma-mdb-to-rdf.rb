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

      meta = JSON.parse(db['meta'])
      name = meta['gemma-wrapper']['meta']['name']
      trait = meta['gemma-wrapper']['meta']['trait']

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
        result.push [chr, pos, snp, af.round(3), se.round(3), effect.round(3), minusLogP.round(2)]
      end
      env.close
      p ["name","trait","chr","pos","af","se","effect","-LogP"]
      if options[:sort]
        $stderr.print("Sorting...\n")
        result = result.sort_by { |rec| rec[rec.size-1] }.reverse
      end

      result.each { |l| p [name,trait] + l }
    end
  end
end # tmpdir
