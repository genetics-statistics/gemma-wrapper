#!/usr/bin/env ruby
#
# Convert a geno file to lmdb. The lmdb file contains 3 tables: (1)
# genotypes as a list of bytes/numbers, (2) markers as a name with
# possibly some metadata -- both indexed on a packed chr+pos. And then
# there is metadata in (3) the info table.
#
# Currently 2 basic storage formats are supported:
#
#    (1) 4-byte floating point and (2) 1-byte translated values.
#
# It is easy to support more formats by adapting --pack and --gpack.
#
# Four translations are supported out of the box:
#
#    Gf   = genotypes as 4-byte floats where missing values are NaN
#    Gb   = bytes that basically map any value between 0.0..2.0 to 0..254. 255 is for missing values
#    G0_1 = similar to Gb but assumes input values 0, 0.5 and 1.0. 255 is for missing values
#    G0_2 = similar to Gb but assumes input values 0, 1 and 2. 255 is for missing values
#
# These are supported by pangemma. It is possible to try other translations, but it may need adaptation
# of pangemma.
#
# Example of doing a G0_2:
#
#   ./bin/geno2mdb.rb BXD.geno.bimbam --eval '{"0"=>0,"1"=>1,"2"=>2,"NA"=>255}' --pack 'C*' --anno snps.txt.mdb --geno-json BXD.geno.json
#   ./bin/geno2mdb.rb BXD.geno.bimbam --eval G0_2 --anno snps.txt.mdb --geno-json BXD.geno.json
#
# The geno-json is optional and basically adds metadata to the DB.
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

BATCH_SIZE = 10_000
CHRPOS_PACK="S>L>L>" # chr, pos, line. L is uint32, S is uint16 - total 64bit

# Translation tables to char/int
Gf = "to_float_or_nan(g)"
Gb = "(g=='NA' ? 255 : (g.to_f*127.0)).to_i"
G0_1 = { "0"=> 0, "0.5"=> 1, "1" => 2, "NA" => 255 }
G0_2 = { "0"=> 0, "1"=> 1, "2" => 2, "NA" => 255 }

Gfmsg =   { type: "Gf", text: "transform to float", geval: Gf, pack: 'f*' }
Gbmsg =   { type: "Gb", text: "transform to byte (255 is NA)", geval: Gb, pack: 'C*' }
G0_1msg = { type: "G0_1", text: "transform 0,0.5,1,NA to byte values 0,1,2,255", eval: G0_1, pack: 'C*' }
G0_2msg = { type: "G0_2", text: "Transform 0,1,2,NA to byte values 0,1,2,255", eval: G0_2, pack: 'C*' }


options = { show_help: false, input: "BIMBAM", geval: "Gf", pack: Gfmsg[:pack] }

opts = OptionParser.new do |o|
  o.banner = "\nUsage: #{File.basename($0)} [options] filename(s)"

  o.on('-i','--input TYPE', ['BIMBAM'], 'input type BIMBAM (default)') do |type|
    options[:input] = type
  end

  o.on('-e','--eval EVAL',String, "eval conversion - note the short cut methods G0_1,G0_2 (default is G0_2)
                                     Example: --eval {\"0\"=>0,\"1\"=>1,\"2\"=>2,\"NA\"=>255} or --eval G0_1

#{Gfmsg}
#{Gbmsg}
#{G0_1msg}
#{G0_2msg}
    ") do |eval|
    options[:format] = eval
    case eval
    when 'Gf'
      options[:geval] = Gf
      options[:pack] = Gfmsg[:pack]
    when 'Gb'
      options[:geval] = Gb
      options[:pack] = Gbmsg[:pack]
    else
      options[:eval] = eval
    end
  end

  o.on('-g','--geval EVAL',String, 'generic eval conversion without assuming it is a hash ([g] is not attached).
                                     Example: --geval {"0"=>0,"1"=>1,"2"=>2,"NA"=>255}[g]
       ') do |geval|
    options[:format] = geval
    options[:geval] = geval
  end

  o.on('-p','--pack PACK',String, 'pack result
                                     Example: "C*"
       ') do |pack|
    options[:pack] = pack
  end

  o.on('-p','--gpack PACK',String, 'generic pack result using your own method
                                     Example: --gpack \'l.pack("C*")\'
         --gpack "l.each_slice(4).map { |slice| slice.map.with_index.sum {|val,i| val << (i*2) } }.pack(\"C*\")"
       ') do |pack|
    options[:gpack] = pack
  end

  o.on('--anno FILEN', 'mdb annotation file') do |anno|
    options[:anno] = anno
    raise "Annotation input file #{anno} does not exist" if !File.exist?(anno)
  end

  o.on('--geno-json JSON', 'JSON gn-geno-to-gemma annotation file') do |json|
    options[:geno_json] = json
  end

  o.on_tail('--order', 'Force rewrite of ordere key-value store') do
    options[:order] = true
  end

  # o.on("-v", "--verbose", "Run verbosely") do |v|
  #   options[:verbose] = true
  # end

  o.on("-d", "--debug", "Show debug messages and keep intermediate output") do |v|
    options[:debug] = true
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
  p options
  exit 1
end

def to_float_or_nan(str)
  Float(str)
rescue ArgumentError
  Float::NAN
end

PACK = if options[:gpack]
         options[:gpack]
       else
         "l.pack(\"#{options[:pack]}\")"
       end
$stderr.print "P_lambda = lambda { |l| #{PACK} }\n"
P_lambda = eval "lambda { |l| #{PACK} }"

EVAL = if options[:geval]
         options[:geval]
       else
         options[:eval]+"[g]"
       end
# G_lambda = eval "lambda { |g| p [EVAL, g, #{EVAL}] ; #{EVAL} }"
G_lambda = eval "lambda { |g| #{EVAL} }"
$stderr.print "G_lambda = lambda { |g| #{EVAL} }\n"

def convert gs, g_func
  res = gs.map { | g | g_func.call(g) }
  # p res
  missing =
    if res[0].is_a?(Float)
      res.count{ |g| g.nan? }
    else
      res.count(255) # byte values use 255 for missing data
    end
  # p [:missing,missing]
  # p res
  # original res.pack(PACK)
  return P_lambda.call(res), missing
end

numsamples =
  if options[:geno_json]
    json = JSON.parse(File.read(options[:geno_json]))
    raise "We need a gemma-geno style JSON!" if not json["type"]=="gn-geno-to-gemma"
    json["numsamples"].to_i
  else
    -1
  end

meta = {
  "type" => "gemma-geno",
  "format" => options[:format],
  "version" => 1.0,
  "eval" => EVAL.to_s,
  "key-format" => CHRPOS_PACK,
  "rec-format" => PACK,
  "geno" => json
}

annofn = options[:anno]
$stderr.print "Reading #{annofn}\n"
marker_env = LMDB.new(annofn, nosubdir: true)
begin
  anno_marker_tab = marker_env.database("marker", create: false)
rescue
  raise "Problem reading annotation file #{annofn}!"
end
keys_ordered = 0
prev_key = ""
cols = -1

ARGV.each_with_index do |fn|
  $stderr.print "Reading #{fn}\n"
  mdb = fn + ".mdb"
  File.delete(mdb) if File.exist?(mdb)
  $stderr.print("Writing lmdb #{mdb}...")
  env = LMDB.new(mdb, nosubdir: true,
                 mapsize: 10**12,
                 maxdbs: 10)
  # maindb      = env.database
  geno        = env.database("geno", create: true)
  geno_marker = env.database("marker", create: true)
  maindb = env.database
  p options

  count = 0
  total_missing = 0
  File.open(fn).each_line.each_slice(BATCH_SIZE) do |batch|
    print "."
    env.transaction() do
      batch.each do |line|
        count += 1
        marker,loc1,loc2,*gs = line.split(/[\s,]+/)
        snpchr = anno_marker_tab[marker]
        raise "Unknown marker #{marker} in #{annofn}!" if !snpchr
        if cols != -1
          raise "Differing amount of genotypes at line #{count}: #{line}" if cols != gs.size
        else
          cols = gs.size
          numsamples = cols if numsamples == -1
          raise "Wrong number of samples in JSON #{numsamples} for #{cols}" if cols != numsamples
        end
        begin
          # key = marker.force_encoding("ASCII-8BIT")
          chr,pos,num = snpchr.unpack(CHRPOS_PACK)
          key = [chr,pos,num].pack(CHRPOS_PACK)
          raise "key error" if not key==snpchr
          keys_ordered += 1 if key >= prev_key
          geno_marker[key] = marker
          fields,missing = # Convert fields to array of values
            case EVAL
            when  "Gf"
              convert(gs, lambda { |g| to_float_or_nan(g) })
            when  "Gb"
              convert(gs, G_lambda)
            when  "G0_1"
              convert(gs, lambda { |g| G0_1[g] })
            when  "G0_2"
              convert(gs, lambda { |g| G0_2[g] })
            else
              convert(gs, G_lambda)
            end
          # p fields
          geno[key] = fields
          prev_key = key
          # track missing data
          total_missing += missing
        rescue TypeError
          raise "Problem at line #{count}: #{line}"
        end
      end
    end
  end
  info = env.database("info", create: true)
  info['numsamples'] = [numsamples].pack("Q") # uint64
  info['nummarkers'] = [geno.size].pack("Q")
  info['meta'] = meta.to_json.to_s
  info['format'] = options[:format].to_s
  info['options'] = options.to_s
  $stderr.print "#{keys_ordered}/#{count} keys are ordered (#{((1.0*keys_ordered/count)*100.0).round(0)}%)\n"
  $stderr.print "#{count-geno.size}/#{geno.size} are duplicate keys!\n"
  $stderr.print "We have #{total_missing} missing values #{(100.0*total_missing/(count*cols)).round(0)}%!\n"
  fn_o = fn + ".order.mdb"
  fn = fn + '.mdb'
  File.delete(fn_o) if File.exist?(fn_o)
  if options[:order] # rewrite ordered store
    $stderr.print "Reordering store #{fn}\n"
    o_env = LMDB.new(fn_o, nosubdir: true, mapsize: 10**9)
    o_geno = o_env.database('geno', create: true)
    o_geno_marker = o_env.database('marker', create: true)
    o_info = o_env.database('info', create: true)
    geno.each do | k,v |
      o_geno[k] = v
    end
    geno_marker.each do | k, v |
      o_geno_marker[k] = v
    end
    info.each do | k, v |
      o_info[k] = v
    end
    o_env.close
  end
  env.close
  File.rename(fn_o,fn) if options[:order]
  $stderr.print "Testing #{fn}\n"
  test_env = LMDB.new(fn, nosubdir: true)
  test_info = test_env.database('info', :create=>false)
  meta2 = test_info.get "meta"
  print meta2,"\n"
  # 1       rs13476251      174792257
  geno_tab = test_env.database('geno', :create=>false)
  marker_tab = test_env.database('marker', :create=>false)
  key = [1,174792257,0].pack(CHRPOS_PACK)
  p marker_tab[key] if marker_tab[key]
end
