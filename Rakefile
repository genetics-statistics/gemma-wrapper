# encoding: utf-8
#
# Run tests with, for example
#
#   env GEMMA_COMMAND=../gemma/bin/gemma rake test

require 'rubygems'
require 'rake'

task default: %w[test]

task :test do
  ruby "bin/gemma-wrapper --json --force -- \
        -g test/data/input/BXD_geno.txt.gz \
        -p test/data/input/BXD_pheno.txt \
        -a test/data/input/BXD_snps.txt \
        -gk \
        -debug > K0.json"
  K0 = File.read("K0.json")
  fail "Wrong Hash in #{K0}" if K0 !~ /1b700de28f242d561fc6769a07d88403764a996f/
  fail "Expected error is 0 in #{K0}" if K0 !~ /errno\":0/
  fail "Test failed" if $? != 0
  ruby "bin/gemma-wrapper --json --input K0.json -- \
        -g test/data/input/BXD_geno.txt.gz \
        -p test/data/input/BXD_pheno.txt \
        -c test/data/input/BXD_covariates2.txt \
        -a test/data/input/BXD_snps.txt \
        -lmm 2 -maf 0.1 \
        -debug > GWA0.json"
  gwa0 = File.read("GWA0.json")
  fail "Wrong Hash in #{gwa0}" if gwa0 !~ /9e411810ad341de6456ce0c6efd4f973356d0bad/
  fail "Expected cache hit in #{gwa0}" if gwa0 !~ /cache_hit\":true/
  fail "Test failed" if $? != 0
  ruby "bin/gemma-wrapper --debug --json --force \
        --loco --chromosomes 1,2,3,4 -- \
        -g test/data/input/BXD_geno.txt.gz \
        -p test/data/input/BXD_pheno.txt \
        -a test/data/input/BXD_snps.txt \
        -gk -debug > KLOCO1.json"
  kloco1 = File.read("KLOCO1.json")
  fail "Wrong Hash in #{kloco1}" if kloco1 !~ /1b700de28f242d561fc6769a07d88403764a996f/
  fail "Expected error is 0 in #{kloco1}" if kloco1 !~ /errno\":0/
  fail "Test failed" if $? != 0
  # run again for cache hits
  ruby "bin/gemma-wrapper  --json  \
        --loco --chromosomes 1,2,3,4 -- \
        -g test/data/input/BXD_geno.txt.gz \
        -p test/data/input/BXD_pheno.txt \
        -a test/data/input/BXD_snps.txt \
        -gk -debug > KLOCO2.json"
  kloco2 = File.read("KLOCO2.json")
  fail "Wrong Hash in #{kloco2}" if kloco2 !~ /1b700de28f242d561fc6769a07d88403764a996f/
  fail "Expected cache hit in #{kloco2}" if kloco2 !~ /cache_hit\":true/
  fail "Test failed" if $? != 0
  ruby "bin/gemma-wrapper --json --force --loco --input KLOCO1.json -- \
        -g test/data/input/BXD_geno.txt.gz \
        -p test/data/input/BXD_pheno.txt \
        -c test/data/input/BXD_covariates2.txt \
        -a test/data/input/BXD_snps.txt \
        -lmm 2 -maf 0.1 \
        -debug > GWA1.json"
  gwa1 = File.read("GWA1.json")
  fail "Wrong Hash in #{gwa1}" if gwa1 !~ /9e411810ad341de6456ce0c6efd4f973356d0bad/
  fail "Test failed" if $? != 0
  # and run again
  ruby "bin/gemma-wrapper --json --loco --input KLOCO2.json -- \
        -g test/data/input/BXD_geno.txt.gz \
        -p test/data/input/BXD_pheno.txt \
        -c test/data/input/BXD_covariates2.txt \
        -a test/data/input/BXD_snps.txt \
        -lmm 2 -maf 0.1 \
        -debug > GWA2.json"
  fail "Test failed" if $? != 0
  gwa2 = File.read("GWA2.json")
  fail "Wrong Hash in #{gwa2}" if gwa2 !~ /9e411810ad341de6456ce0c6efd4f973356d0bad/
  fail "Expected cache hit in #{gwa2}" if gwa2 !~ /cache_hit\":true/
end

require 'rdoc/task'
Rake::RDocTask.new do |rdoc|
  version = File.exist?('VERSION') ? File.read('VERSION') : ""

  rdoc.rdoc_dir = 'rdoc'
  rdoc.title = "bio-gemma-wrapper #{version}"
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/**/*.rb')
end
