#!/usr/bin/env ruby
#
# Convert a marker (snp) annotation file to lmdb - ends up being
# larger, but maybe faster.  The key is stored as a packed "S>L>L>"
# for [chr,pos,line] where chr is mapped to a number 0..20 and maps
# X,Y,M to its ASCII values. line ascertains the key is unique. The
# data is simply the marker name as a string.
#
# If you get a compatibility error in guix you may need an older
# Ruby. Otherwise you can do:
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

BATCH_SIZE=10_000 # increasing does not really speed things up
CHRPOS_PACK="S>L>L>" # L is uint32, S is uint16 - total 64bit

options = { show_help: false }

opts = OptionParser.new do |o|
  o.banner = "\nUsage: #{File.basename($0)} [options] filename(s)"

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

dup_count = 0
ARGV.each do |fn|
  $stderr.print "Reading #{fn}..."
  mdb = fn + ".mdb"
  File.delete(mdb) if File.exist?(mdb)
  $stderr.print("Writing lmdb #{mdb}...")
  env = LMDB.new(mdb, nosubdir: true, nosync: true, mapsize: 10**9)
  maindb = env.database
  marker_tab = env.database("marker", create: true) # store chrpos->marker
  info = env.database("info", create: true)
  meta = {
    "type" => "gemma-anno",
    "version" => 1.0,
    "key-format" => CHRPOS_PACK,
    "rec-format" => "string",
  }
  info["meta"] = meta.to_json.to_s

  count = 0
  File.open(fn).each_line.each_slice(BATCH_SIZE) do |batch|
    print "."
    env.transaction() do
      batch.each do |line|
        count += 1
        name,pos,chr = line.split(/[\s,]+/)
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
        chrposdup = [chr_c,pos_i,count].pack(CHRPOS_PACK) # count handles duplicates
        marker_tab[chrposdup] = name
      end
    end
  end
  $stderr.print "\n#{marker_tab.size}/#{count} records written\n"
  env.close
end
