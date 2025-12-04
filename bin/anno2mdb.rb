#!/usr/bin/env ruby
#
# Convert a (snp) annotation file to lmdb - ends up being larger, but maybe faster
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

CHRPOS_PACK="L>L>"

options = { show_help: false }

opts = OptionParser.new do |o|
  o.banner = "\nUsage: #{File.basename($0)} [options] filename(s)"

  # o.on('-a','--anno FILEN', 'Annotation file') do |anno|
  #   options[:anno] = anno
  #   raise "Annotation input file #{anno} does not exist" if !File.exist?(anno)
  # end

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

X='X'.ord
Y='Y'.ord
M='M'.ord

ARGV.each do |fn|
  $stderr.print "Reading #{fn}\n"
  mdb = fn + ".mdb"
  File.delete(mdb) if File.exist?(mdb)
  $stderr.print("Writing lmdb #{mdb}...\n")
  env = LMDB.new(mdb, nosubdir: true, mapsize: 10**9)
  maindb = env.database
  chrpos_tab = env.database("chrpos", create: true, integerkey: true, dupsort: true)
  marker_tab = env.database("marker", create: true) # store reversed marker -> chrpos

  count = 0
  File.open(fn).each_line do |line|
    count += 1
    snp,pos,chr = line.split(/[\s,]+/)
    location = "#{chr}:#{pos}"
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
    begin
      pos_i = Integer(pos)
    rescue ArgumentError, TypeError
      pos_i = 0 # set anything unknown to position zero
    end
    chrpos = [chr_c,pos_i].pack(CHRPOS_PACK)
    chrpos_tab[chrpos] = snp
    marker_tab[snp] = chrpos
  end
  $stderr.print "#{chrpos_tab.size}/#{count} records written\n"
  env.close

end
