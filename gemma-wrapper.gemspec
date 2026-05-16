Gem::Specification.new do |s|
  s.name        = 'bio-gemma-wrapper'
  s.version     = File.read('VERSION')
  s.summary     = "GEMMA with LOCO and permutations"
  s.description = "GEMMA wrapper adds LOCO and permutation support. Also runs in parallel and caches K between runs with LOCO support"
  s.authors     = ["Pjotr Prins"]
  s.email       = 'pjotr.public01@thebird.nl'
  s.files       = Dir['bin/*'].reject { |f| File.directory?(f) } +
                  Dir['lib/**/*.rb'] +
                  Dir['test/**/*'].reject { |f| File.directory?(f) } +
                  ["Rakefile",
                   "Gemfile",
                   "LICENSE.txt",
                   "README.md",
                   "RELEASE_NOTES.md",
                   "VERSION",
                   "gemma-wrapper.gemspec"]
  s.executables = Dir['bin/*'].reject { |f| File.directory?(f) }.map { |f| File.basename(f) }
  s.homepage    =
    'https://github.com/genetics-statistics/gemma-wrapper'
  s.license       = 'GPL3'
  s.require_paths = ["lib"]
  s.required_ruby_version = Gem::Requirement.new(">= 2.0.0")

  # bin/gemma-wrapper itself only uses the Ruby stdlib (json, optparse,
  # tmpdir, socket, etc.) plus lib/lock.rb which ships in the gem.  The
  # helper scripts under bin/ pull in extra gems; declare them as
  # optional runtime deps so `gem install bio-gemma-wrapper` produces a
  # working install for the *mdb* / RDF helpers as well.  The lmdb gem
  # backs anno-mdb-to-rdf, anno2mdb, gemma-mdb-to-rdf and geno2mdb;
  # rdf + rdf-vocab back the *-to-rdf helpers.
  s.add_runtime_dependency 'lmdb', '~> 0.6'
  s.add_runtime_dependency 'rdf', '~> 3.0'
  s.add_runtime_dependency 'rdf-vocab', '~> 3.0'

  s.add_development_dependency 'rake', '~> 13.0'
end
