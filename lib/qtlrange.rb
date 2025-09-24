# This module handles a list of chromosomes and QTL ranges, tracking SNPs

MAX_SNP_DISTANCE_BPS = 50_000_000.0
MAX_SNP_DISTANCE = MAX_SNP_DISTANCE_BPS/10**6
LOD_MIN = 4.0

module QTL

  class QLocus
    attr :id,:chr,:pos,:af,:lod
    def initialize id,chr,pos,af,lod
      @id = id
      @chr = chr
      @pos = pos
      @af = af
      @lod = lod
    end
  end

  # Track one QTL using a range
  class QRange
    attr_reader :chr,:min,:max,:snps,:max_af,:lod
    def initialize locus
      @snps = [locus.id]
      @chr = locus.chr
      @min = locus.pos
      @max = locus.pos
      @max_af = locus.af
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
      @max_af=locus.af if @max_af==nil or locus.af > @max_af
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
      return true if qtl.min >= @min and qtl.max <= @max # qtl falls with bounds
      return true if qtl.min <= @min and qtl.max >= @min # qtl over left boundary
      return true if qtl.min <= @max and qtl.max >= @max # qtl over right boundary
      false
    end

    def inspect
      if @max_af != nil
        "#<QRange Chr#{@chr} 𝚺#{snps.size} #{self.min}..#{self.max} AF=#{@max_af} LOD=#{@lod}>"
      else
        "#<QRange Chr#{@chr} 𝚺#{snps.size} #{self.min}..#{self.max} LOD=#{@lod}>"
      end
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
"""
          print "    gnt:qtlAF       #{qtl.max_af} ;\n" if qtl.max_af
          print "    gnt:qtlLOD      #{qtl.lod.max} .\n"
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
        if print_rdf # Output QTL to RDF
          qtlid = gnqtlid(id,qtl2)
          if qtl2.lod.max >= LOD_MIN
            print "#{qtlid} a gnt:newQTL .\n"
          end
        end
      end
      $stderr.print ["#{name}: NO #{set1.method} results, new QTL(s) #{set2.method} Chr #{chr}!",qtls2],"\n"
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
            if qtl2.lod.max >= LOD_MIN
              print "#{qtlid} a gnt:newQTL .\n"
            end
          end
          $stderr.print ["#{name}: NO #{set1.method} match, QTL #{set2.method} Chr #{chr}!",qtl2],"\n"
        end
      end
    end
  end
end
