#spotNum,x0,y0,pixels,sum,maxPixels,pixelsHit,pixelsMissed,fractColoc,tag
#0,1258.04,20.5926,27,14618,652,0,27,0, accepted

class CellResult
   @name
   @n
   @pixels_sum
   @pixels_ave
   @tot_sum
   @tot_ave
   @n_coloc
   @n_missed
   @pixels_coloc_sum
   @pixels_coloc_ave
   @pixels_missed_sum
   @pixels_missed_ave
   @tot_coloc_sum
   @tot_coloc_ave
   @tot_missed_sum
   @tot_missed_ave
   @coloc_thresh
   
   def calc()
      @pixels_ave        = @pixels_sum / @n
      @tot_ave           = @tot_sum / @n
      @pixels_coloc_ave  = @pixels_coloc_sum / @n_coloc
      @pixels_missed_ave = @pixels_missed_sum / @n_missed
      @tot_coloc_ave     = @tot_coloc_sum / @tot_coloc
      @tot_missed_ave    = @tot_missed_sum / @tot_missed
   end

   def initialize(name, coloc_thresh)
      @name = name
      @n = @pixels_sum = @pixels_ave = @tot_sum = @tot_ave = 0.0
      @n_coloc = @n_missed = 0.0
      @pixels_coloc_sum = @pixels_missed_sum = 0.0
      @tot_coloc_sum = @tot_missed_sum = 0.0
      @tot_coloc = @tot_missed = 0.0
      @coloc_thresh = coloc_thresh
   end

   def addLine(line)
      return if line.length < 10
      @n += 1
      @pixels_sum += line[3].to_f
      @tot_sum    += line[4].to_f
      if @coloc_thresh <  line[8].to_f then
         @n_coloc += 1
         @pixels_coloc_sum += line[3].to_f
         @tot_coloc_sum    += line[4].to_f
      else
         @n_missed += 1
         @pixels_missed_sum += line[3].to_f
         @tot_missed_sum    += line[4].to_f
      end
   end

   def report_header() 
      return "Name,TotalSpots,NumColoc,NumMissed,MassColoc,MassMissed"
   end

   def report_line()
      self.calc()
      "#{@name},#{@n},#{@n_coloc},#{@n_missed},#{@tot_coloc_sum},#{@tot_missed_sum}"
   end
end

class LineEq
   @m
   @b
   @x1
   @x2
   @y1
   @y2

   def initialize(x1, y1, x2, y2)
       @x1 = x1.to_f
       @y1 = y1.to_f
       @x2 = x2.to_f
       @y2 = y2.to_f
       @m = (@y2 - @y1) / (@x2 - @x1)
       @b = @y1 - @m*x1
   end

   def y(x)       # not defined for horizontal line
     if @y1 = @y2 then
        return y1        # line horizontal, y never varies
     end
     @m*x + @b
   end
 
   def x(y)       # not defined for vertical line
     if @x1 = @x2 then
        return @x1        # line vertical, x never varies
     end
     (y - @b) / @m
   end

end

class Edgemap
   @polygon        # list of line segments
   @edges
   @ymax

   def printEdgemap()
      puts "Number of edges = #{@edges.length}"
      i=0
      puts " ----- STARTING -----"
      @edges.each { | e |
         if !e || (e==[]) then
           puts "[nil]"
         else
           puts "#{i}:#{e}"
         end
         i += 1
      }
      puts " ----- ENDING -----"
   end

   def addedges(start, finish) 
      x0, y0 = start
      x1, y1 = finish
      if y0 == y1 then
         #@edges[y0]   = [] if !@edges[y0]
         @edges[y0] << x0 << x1
         return
      end
      
      # calculate line
      line = LineEq.new(x0,y0,x1,y1)
      if y0 < y1 then
         for i in y0..(y1-1) do
            #puts "#{i}:*"
            #@edges[i] = [] if !@edges[i]
            @edges[i] << line.x(i).round  
         end
      else   # y0 > y1
         for i in (y1-1)..y0 do
            #puts "#{i}:-"
            #@edges[i] = [] if !@edges[i]
            @edges[i] << line.x(i).round  
         end
      end      
   end

   # each line segment is taken as closed-open

   def makeEdgemap(poly, ymax)
      lastpoint = nil
      @edges    = Array.new(ymax,0)
      puts "Length of @edges is #{@edges.length}"
      for i in 0..@edges.length do 
         @edges[i] = []
      end 
      poly.each { | point | 
         point[1] = point[1].to_i
         point[2] = point[2].to_i
         if lastpoint then
            addedges(lastpoint, point)
         end
         lastpoint = point
      }
   end

   def inEdge?(ptlist, x) 
      answer = false
      #puts "#{x} being tested in #{ptlist}"
      return false if !ptlist || (ptlist.length == 0)
      ptlist.each { | p |
         if p < x then
            answer = answer ? false : true
         else
            return answer
         end
      }
      return answer
   end

   def inPoly?(x,y) 
      if @edges[y.to_i]==[] then
         return false
      end
      level = @edges[y.to_i]
      self.inEdge?(level,x.to_i)
   end

   def initialize(polygon, ysize) 
      @polygon = polygon      # not sure why I need to keep it
      # make sure polygon is closed
      if (@polygon[0][0] != @polygon[-1][0]) && (@polygon[0][1] != @polygon[-1][1]) then
          @polygon << @polygon[-1]
      end
      makeEdgemap(@polygon, ysize)
   end

end

#--------------------

def readInCells(filename) 
   cellname = ""
   pointlist= []
   cells    = []
   inFile   = File.new(filename)
   inFile.each_line("\n") { | line |
      l2 = line.split(" ")
      if l2.length == 3 then
         cell = l2[0]
         x    = l2[1].to_i
         y    = l2[2].to_i
         if cellname != cell then
            cells << pointlist if pointlist.length > 0
            pointlist = []
            cellname = cell
         end
         pointlist << [x,y]
      end
   }
   inFile.close
   cells << pointlist if pointlist.length > 0 
end

def readInData(filename) 
   inFile   = File.new(filename)
   spotlist = []
   inFile.each_line("\n") { | line |
      if line =~ /pixels/      # skip header 
         next
      end
      items = line.split(",")
      if items.length >= 10 then
         spotlist << items     # note that many strings need to be converted to int or float
      end
   }
   spotlist
end

# ------ main ------
if ARGV.length != 3 then
   puts "coloc.rb tip47a.txt tip47af.dat tip47a.results"
   return
end
#cellFile    = "tip47a.txt"
cellFile    = ARGV[0]
cells       = readInCells(cellFile)
datFile     = ARGV[1]
#datFile     = "tip47af.dat"
dat         = readInData(datFile)
#outfilename = "out_test.dat"
outfilename = ARGV[2]
#if File.exists(outfilename) then
   outfile = File.new(outfilename,"a")
#else
#   outfile = File.new(outfilename,"w")
#end
headerWritten = false
cellResults = []
counter     = 0
# -- now divide data between cells ignoring points not in cells
cells.each { | cell |
  name = "#{cellFile}:cell#{counter.to_s}"
  counter += 1
  cellResult = CellResult.new(name, 0.3)
  if !headerWritten then
      outfile << cellResult.report_header() << "\n" 
      headerWritten = true
  end
  edges      = Edgemap.new(cell, 2048) 
  #edges.printEdgemap
  dat.each { | d | 
     if edges.inPoly?(d[1].to_f.round, d[2].to_f.round) then
       cellResult.addLine(d)
     else
       next
       #puts "#{d} failed."
     end
  }
  outfile << cellResult.report_line() << "\n"
}
outfile.close()
