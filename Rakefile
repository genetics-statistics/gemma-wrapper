# encoding: utf-8
#
# Run tests with, for example
#
#   env GEMMA_COMMAND=../gemma/bin/gemma rake test

require 'rubygems'
require 'rake'

task default: %w[test]

task :test do
  ruby "bin/gemma-wrapper  --json --force \
        --loco 1,2,3,4 -- \
        -g test/data/input/BXD_geno.txt.gz \
        -p test/data/input/BXD_pheno.txt \
        -a test/data/input/BXD_snps.txt \
        -gk -debug > K.json"
  fail "Test failed" if $? != 0
  # run again for cache hits
  ruby "bin/gemma-wrapper  --json  \
        --loco 1,2,3,4 -- \
        -g test/data/input/BXD_geno.txt.gz \
        -p test/data/input/BXD_pheno.txt \
        -a test/data/input/BXD_snps.txt \
        -gk -debug > K2.json"
  fail "Test failed" if $? != 0
  ruby "bin/gemma-wrapper --json --force --loco --input K.json -- \
        -g test/data/input/BXD_geno.txt.gz \
        -p test/data/input/BXD_pheno.txt \
        -c test/data/input/BXD_covariates2.txt \
        -a test/data/input/BXD_snps.txt \
        -lmm 2 -maf 0.1 \
        -debug > GWA.json"
  fail "Test failed" if $? != 0
  # and run again
  ruby "bin/gemma-wrapper --json --loco --input K.json -- \
        -g test/data/input/BXD_geno.txt.gz \
        -p test/data/input/BXD_pheno.txt \
        -c test/data/input/BXD_covariates2.txt \
        -a test/data/input/BXD_snps.txt \
        -lmm 2 -maf 0.1 \
        -debug > GWA2.json"
  fail "Test failed" if $? != 0
end

require 'rdoc/task'
Rake::RDocTask.new do |rdoc|
  version = File.exist?('VERSION') ? File.read('VERSION') : ""

  rdoc.rdoc_dir = 'rdoc'
  rdoc.title = "bio-gemma-wrapper #{version}"
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/**/*.rb')
end
