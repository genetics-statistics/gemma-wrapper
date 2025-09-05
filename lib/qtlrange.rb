# This module handles a list of chromosomes and QTL ranges, tracking SNPs

=begin

          r = Range.new(1.1,2.0)
          irb(main):013:0> [4.0,r.max].max
          => 4.0
          irb(main):014:0> [1.0,r.max].max
          => 2.0

=end


module QTL


  # Track one QTL using a range
  class QRange < Range
    attr_reader :chr,:snps
    def initialize id, chr, pos
      @snps = [id]
      @chr = chr
      super(pos,pos)
    end

    def inspect
      "#<QRange 𝚺#{snps.size} #{self.min}..#{self.max}>"
    end
  end

  # Track all ranges
  class QRanges
    attr_reader :chromosome
    def initialize
      @chromosome = {}
      @snps = {}
    end

    def add_snp snp_id, chr, pos
      @chromosome[chr] = [] if not @chromosome.has_key? chr
      # qtl = QRange(snp_id,snp)
      ranges = @chromosome[chr]
      hit = QRange.new(snp_id,chr,pos)
      covered = false
      ranges.each do |range|
        if range.include?(pos)
          covered = true
        end
      end
      ranges.append(hit) if not covered

        # if not ranges.has_key?(chr)
        #   ranges[chr] = [ Range.new(pos,pos) ]
        # else
        #   covered = false
        #   ranges.values.each do | range |
        #     p range
        #     covered = true  if range.include?(pos)
        #   end
        #   ranges[chr].append Range.new(pos,pos) if not covered
        # end
    end
  end

end
