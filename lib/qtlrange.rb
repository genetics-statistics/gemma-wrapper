# This module handles a list of chromosomes and QTL ranges, tracking SNPs

MAX_SNP_DISTANCE_BPS = 50_000_000.0
MAX_SNP_DISTANCE = MAX_SNP_DISTANCE_BPS/10**6

module QTL

  class QLocus
    attr :id,:chr,:pos,:lod
    def initialize id,chr,pos,lod=nil
      @id = id
      @chr = chr
      @pos = pos
      @lod = lod
    end
  end

  # Track one QTL using a range
  class QRange
    attr_reader :chr,:min,:max,:snps,:lod
    def initialize locus
      @snps = [locus.id]
      @chr = locus.chr
      @min = locus.pos
      @max = locus.pos
      @lod = Range.new(locus.lod,locus.lod)
    end

    def add locus
      chr = locus.chr
      raise "Chr #{chr} mismatched for range #{self}" if chr != @chr
      @snps.append(locus.id)
      @min = [locus.pos,@min].min
      @max = [locus.pos,@max].max
      test_lod = locus.lod
      return if test_lod == nil
      if @lod.min == nil or @lod.max == nil
        @lod = Range.new(test_lod,test_lod)
      else
        @lod = Range.new([test_lod,@lod.min].min,[test_lod,@lod.max].max)
      end
    end

    def in_range? locus
      pos = locus.pos
      pos > @min - MAX_SNP_DISTANCE and pos < @max + MAX_SNP_DISTANCE
    end

    # See if a QRange overlaps with one of ours
    def overlaps? qtl
      return true if qtl.min > @min and qtl.max < @max
      return true if qtl.min < @min and qtl.max > @min
      return true if qtl.min < @max and qtl.max > @max
      false
    end

    def inspect
      "#<QRange Chr#{@chr} 𝚺#{snps.size} #{self.min}..#{self.max} LOD=#{@lod}>"
    end
  end

  # Track all ranges
  class QRanges
    attr_reader :chromosome,:name,:method
    def initialize name, method=""
      @chromosome = {}
      @name = name
      @method = method
    end

    def add_locus locus
      chr = locus.chr
      @chromosome[chr] = [] if not @chromosome.has_key? chr
      ranges = @chromosome[chr]
      covered = false
      nrange = QRange.new(locus)
      ranges.each do |range|
        if range.in_range?(locus)
          range.add(locus)
          covered = true
        end
      end
      ranges.append(nrange) if not covered
      # make sure they are ordered
      @chromosome[chr] = ranges.sort_by { |r| r.min }
    end

    def print_rdf(id)
      chromosome.sort.each do |chr,r|
        r.each do | qtl |
          qtlid = gnqtlid(id,qtl)
          print """
#{qtlid}
    gnt:mappedQTL   #{id};
    rdfs:label      \"GEMMA BXDPublish QTL\";
    gnt:qtlChr      \"#{chr}\";
    gnt:qtlStart    #{qtl.min} ;
    gnt:qtlStop     #{qtl.max} ;
    gnt:qtlLOD      #{qtl.lod.max} .
"""
          qtl.snps.each do |snp|
            print """#{qtlid} gnt:mappedSnp #{gnid(snp)} .
"""
          end
        end
      end
    end

    def to_s
      "[#{@name},#{@method}] =>{" + chromosome.sort.map{|k,v| "#{k.inspect}=>#{v.inspect}"}.join(", ") + "}"
    end

  end
end

# Take two QRanges and say what is happening
def qtl_diff(id, set1, set2, print_rdf)
  $stderr.print set1,"\n"
  $stderr.print set2,"\n"
  name = set1.name
  # we check by chromosome
  set2.chromosome.each do | chr,qtls2 |
    qtls1 = set1.chromosome[chr]
    if qtls1 == nil
      qtls2.each do | qtl2 |
        if print_rdf
          qtlid = gnqtlid(id,qtl2)
          print "#{qtlid} a gnt:newlyDiscoveredQTL .\n"
        end
      end
      $stderr.print ["#{name}: NO HK #{set1.method} match for #{set2.method} Chr #{chr} QTL!",qtls2],"\n"
    else
      qtls2.each do | qtl2 |
        match = false
        qtls1.each do | qtl1 |
          match = true if qtl2.overlaps?(qtl1)
          break
        end
        if not match
          if print_rdf
            qtlid = gnqtlid(id,qtl2)
            print "#{qtlid} a gnt:newlyDiscoveredQTL .\n"
          end
          $stderr.print ["#{name}: NO #{set1.method} match for #{set2.method} Chr #{chr} QTL!",qtl2],"\n"
        end
      end
    end
  end
end
