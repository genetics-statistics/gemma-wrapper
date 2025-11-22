#!/usr/bin/env ruby
#
# Convert a geno file to lmdb. Example:
#
#   ./bin/geno2mdb.rb BXD.geno.bimbam --eval '{"0"=>0,"1"=>1,"2"=>2,"NA"=>-1}' --pack 'C*' --geno-json BXD.geno.json
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

options = { show_help: false, input: "BIMBAM", eval: "G0-2", pack: "C*" }

opts = OptionParser.new do |o|
  o.banner = "\nUsage: #{File.basename($0)} [options] filename(s)"

  o.on('-i','--input TYPE', ['BIMBAM'], 'input type BIMBAM (default)') do |type|
    options[:input] = type
  end

  o.on('-e','--eval EVAL',String, 'eval conversion - note the short cut methods G0-1,G0-2 are faster (default is G0-2)') do |eval|
    options[:eval] = eval
  end

  o.on('-p','--pack PACK',String, 'pack result') do |pack|
    options[:pack] = pack
  end

  o.on('--geno-json JSON', 'JSON gn-geno-to-gemma annotation file') do |json|
    options[:geno_json] = json
  end

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
  p options
  exit 1
end

PACK=options[:pack]
# Translation tables to char/int
G0_1 = { "0"=> 0, "0.5"=> 1, "1" => 2, "NA" => 255 }
G0_2 = { "0"=> 0, "1"=> 1, "2" => 2, "NA" => 255 }

G_lambda = eval "lambda { |g| #{options[:eval]} }"

def convert gs, func
  res = gs.map { | g | func.call(g) }
  res.pack(PACK)
end

json = JSON.parse(File.read(options[:geno_json]))

meta = {
  "type" => "gemma-geno",
  "version" => 1.0,
  "eval" => options[:eval].to_s,
  "key-format" => "string",
  "rec-format" => PACK,
  "geno" => json
}



cols = -1
ARGV.each_with_index do |fn|
  $stderr.print "Reading #{fn}\n"
  mdb = fn + ".mdb"
  $stderr.print("lmdb #{mdb}...\n")
  env = LMDB.new(mdb, nosubdir: true, mapsize: 10**9)
  maindb = env.database
  db = env.database(File.basename(mdb), create: true)

  count = 0
  File.open(fn).each_line do |line|
    count += 1
    marker,loc1,loc2,*rest = line.split(/[\s,]+/)
    if cols != -1
      raise "Varying amount of genotypes at line #{count}: #{line}" if cols != rest.size
    else
      cols = rest.size
    end
    begin
      db[marker] =
        case options[:eval]
        when  "G0-1"
          convert(rest, lambda { |g| G0_1[g] })
        when  "G0-2"
          convert(rest, lambda { |g| G0_2[g] })
        else
          convert(rest, G_lambda)
        end
    rescue TypeError
      raise "Problem at line #{count}: #{line}"
    end
  end
  db['meta'] = meta.to_json
  env.close
end

raise "Empty set" if cols == -1

# meta["geno"]["cols"] = cols
print meta.to_json,"\n"
