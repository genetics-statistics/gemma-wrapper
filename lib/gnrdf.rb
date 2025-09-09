def gnid(id)
  id.to_s.sub(/^http:\/\/genenetwork.org\/id\//,"gn:")
end

def gnqtlid(id,qtl)
  "#{id}_QTL_Chr#{qtl.chr}_#{qtl.min.round}_#{qtl.max.round}"
end
