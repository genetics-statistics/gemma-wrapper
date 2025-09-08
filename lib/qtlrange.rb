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

    def inspect
      "#<QRange 𝚺#{snps.size} #{self.min}..#{self.max} LOD=#{@lod}>"
    end
  end

  # Track all ranges
  class QRanges
    attr_reader :chromosome
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

    def inspect
      "[#{@name},#{@method}] =>{" + chromosome.sort.map{|k,v| "#{k.inspect}=>#{v.inspect}"}.join(", ") + "}"
    end
  end
end
