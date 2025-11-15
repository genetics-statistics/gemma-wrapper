# This module handles a list of chromosomes and QTL ranges, tracking SNPs
#
# Basically we have a QTL locus (QLocus) that tracks chr,pos,af and lod for each hit.
# QRange is a set of QLocus which also tracks some stats chr,min,max,snps,max_af,lod.
# It can compute whether two QTL (QRange) overlap.
# Next we have a container that tracks the QTL (QRanges) on a chromosome.
#
# Finally there is a diff function that can show the differences on a chromosome (QRanges) for two mapped traits.

MAX_SNP_DISTANCE_BPS = 50_000_000.0
MAX_SNP_DISTANCE = MAX_SNP_DISTANCE_BPS/10**6
LOD_MIN = 4.0

class Hash
  def hmap(&block)
    Hash[self.map {|k, v| block.call(k,v) }]
  end
end

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
      mirror_af = 1.0-locus.af
      @max_af=mirror_af if @max_af==nil or mirror_af > @max_af
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
      return true if qtl.min >= @min and qtl.max <= @max # qtl falls within bounds    {  <   >  }
      return true if @min >= qtl.min and @max <= qtl.max # qtl falls within bounds    <  {   }  >

      return true if qtl.min <= @min and qtl.max >= @min # qtl over left boundary     {  <   }  >
      return true if @min <= qtl.min and @min >= qtl.max # qtl over left boundary     <  {   >  }
      # return true if qtl.min <= @max and qtl.max >= @max # qtl over right boundary  <  {   >  } dup
      # return true if @max <= qtl.min and @max >= qtl.max # qtl over right boundary  {  <   }  > dup
      false
    end

    def combine! qtl
      @min = [@min,qtl.min].min
      @max = [@max,qtl.max].max
      @max_af = [@max_af,qtl.max_af].max
      @lod = Range.new([qtl.lod.min,@lod.min].min,[qtl.lod.max,@lod.max].max)
      @snps = @snps.append(qtl.snps).flatten
      self
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
      # --- track Qranges
      chr = locus.chr
      @chromosome[chr] = [] if not @chromosome.has_key? chr
      ranges = @chromosome[chr]
      # --- here we look for overlap
      covered = false
      nrange = QRange.new(locus)
      ranges.each do |range|
        if range.in_range?(locus)
          range.add(locus)
          covered = true
          break
        end
      end
      ranges.append(nrange) if not covered
      # make sure they are ordered
      @chromosome[chr] = ranges.sort_by { |r| r.min }
    end

    def pangenome_filter
      @chromosome = @chromosome.hmap{ |chr,qtls|
        [ chr, qtls.delete_if { |qtl| qtl.lod.max < 6.0 or (qtl.lod.max < 7.0 - qtl.snps.size/2)} ]
      }.delete_if { |chr, qtls| qtls.empty? }
    end

    def rebin
      # because we grow bins there may still be some overlap
      new_chromosome = {}
      @chromosome.each do |chr, qtls|
        # make sure the bins are sorted
        prev = nil
        new_qtls = []
        qtls.each do | qtl|
          if prev and qtl.overlaps? prev
            new_qtls[-1] = qtl.combine!(prev)
          else
            new_qtls.push qtl
          end
          prev = qtl
        end
        new_chromosome[chr] = new_qtls
      end
      @chromosome = new_chromosome
    end

    def gnqtlid(id,qtl)
      "#{id}_QTL_Chr#{qtl.chr}_#{qtl.min.round}_#{qtl.max.round}"
    end

    def print_rdf(rdf_trait_id)
      chromosome.sort.each do |chr,r|
        r.each do | qtl |
          tid = gnid(rdf_trait_id)
          qtl_id = gnqtlid(tid,qtl)

          print """
#{qtl_id}
    gnt:mappedQTL   #{tid};
    rdfs:label      \"GEMMA BXDPublish QTL\";
    gnt:qtlChr      \"#{chr}\";
    gnt:qtlStart    #{qtl.min} ;
    gnt:qtlStop     #{qtl.max} ;
"""
          print "    gnt:qtlAF       #{qtl.max_af} ;\n" if qtl.max_af
          print "    gnt:qtlLOD      #{qtl.lod.max} .\n"
          qtl.snps.each do |snp|
            print """#{qtl_id} gnt:mappedSnp #{gnid(snp)} .
"""
          end
        end
      end
    end

    def to_s
      # p chromosome
      "[#{@name},#{@method}] =>{" + chromosome.sort_by { |r| r[0].rjust(3) }.map{|k,v| "\n#{k.inspect}=>#{v.inspect}"}.join(", ") + "}"
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
