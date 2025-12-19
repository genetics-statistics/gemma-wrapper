#!/usr/bin/env ruby
#
# Read anno mdb (created by anno2mdb.rb) and write to RDF.
#
# Pjotr Prins (c) 2025

require 'tmpdir'
require 'json'
require 'optparse'

options = { show_help: false }

opts = OptionParser.new do |o|
  o.banner = "\nUsage: #{File.basename($0)} [options] filename"

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

def rdf_normalize(uri)
  uri.gsub(/[-<+>]/,"_")
end

ARGV.each do |traitid|
  $stderr.print "Reading #{traitid}\n"
  cmd = "curl http://127.0.0.1:8091/dataset/bxd-publish/values/#{traitid}.json"
  $stderr.print cmd
  buf = `#{cmd}`
  samples = JSON.parse(buf)
  # p json
  samples.each do | k,v |
    sample = rdf_normalize(k.capitalize)
  print """gn:traitBxd_#{traitid} gnt:sample_id gn:#{sample} .
"""
  end
end
