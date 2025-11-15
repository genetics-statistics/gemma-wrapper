def gnid(id)
  id.to_s.sub(/^http:\/\/genenetwork.org\/id\//,"gn:")
end
