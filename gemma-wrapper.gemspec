Gem::Specification.new do |s|
  s.name        = 'bio-gemma-wrapper'
  s.version     = File.read('VERSION')
  s.summary     = "GEMMA with LOCO and permutations"
  s.description = "GEMMA wrapper adds LOCO and permutation support. Also runs in parallel and caches K between runs with LOCO support"
  s.authors     = ["Pjotr Prins"]
  s.email       = 'pjotr.public01@thebird.nl'
  s.files       = ["bin/gemma-wrapper",
                   "lib/lock.rb",
                   "Rakefile",
                   "Gemfile",
                   "LICENSE.txt",
                   "README.md",
                   "VERSION",
                   "gemma-wrapper.gemspec"
                  ]
  s.executables = ['gemma-wrapper']
  s.homepage    =
    'https://github.com/genetics-statistics/gemma-wrapper'
  s.license       = 'GPL3'
  s.require_paths = ["lib"]
  s.required_ruby_version = Gem::Requirement.new(">= 2.0.0")
end
