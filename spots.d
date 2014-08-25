import std.stdio;
import std.string;
import image;

/*
 *
 * To do:  generic exceptions thrown for now.  Make more specific.
 *
 */

/* This file contains the spots, kspots, nspots, coloc and phSpots routines + their supporting subroutines 
 *     -- spots isolates spots as contigous areas above threshold, or fraction of max int produces report file
 *     -- kspots is famous Fred Maxfield program that gradually raises fraction to final value
 *     -- coloc sees if spots overlap other spots in reference image, deletes otherwise
 *     -- pHspots produces ratios for pH measurements after applying kspots and deleting saturated spots
 *
 *  *** Routines to implement are  coloc, pHspots, edgeSpots
 *  *** A related routine will be threeDspots
 */


struct Point {
  int x;
  int y;
 }


/* initializeQueue(), enq() and deq() provide a queue as a ring buffer for spots
 * the basic algorithm is mark and add pixel, enqueue the locations of all its connected
 * neighbors, then dequeue new pixel to add till the queue is empty.  This is done iteratively
 * to enforce an outward expanding front, so much less intermediate storage is needed than 
 * with the recursive "depth-first" alternative implementation.
 */

struct SpotQueue{
  enum MAXQUEUE = 20000;
  Point spotQueue[MAXQUEUE];
  size_t head = 0;
  size_t tail = 0;
  bool full = false;
  bool empty = true;
  
  void reset() {   // allow reuse of an old queue
     head  = tail = 0;
     full  = false;
     empty = true;
  }
  void enq(Point pixel) {
    if (full) {
       writefln("Head = %d.  Tail = %d.", head, tail);
       throw new Exception("Spot queue full.");
       // later, expanding the array will be an option, too lazy now
    }
    head += 1;
    head = head % MAXQUEUE;
    spotQueue[head] = pixel;
    if ((head == tail)&&(!empty)) {
      full = true;
      writefln("Head = %d.  Tail = %d.", head, tail);
      throw new Exception("Spot queue full.");
    }
    empty = false;
    return;
  }

  Point deq() {
     Point temp;
     if (empty) {
        // If this happens, it is a programming error on part of the calling routine
        throw new Exception("Spot queue empty.");
     }
     tail += 1;
     tail = tail % MAXQUEUE;
     if (tail == head) {
         empty = true;
     } else {
         empty = false;
     }
     return spotQueue[tail];
  }

  bool isEmpty() { return empty; }
} 

struct SpotInfo {
  int spotNo = 0;
  // negative values below  aren't calculated until spot is complete
  double x    = -1;    // centroid
  double y    = -1;
  int x1      = -1;    // bounding box that localizes all pixels
  int y1      = -1;
  int x2      = -1;
  int y2      = -1;
  int sum     = 0;
  int pixels  = 0;
  int sumxs   = 0;
  int sumys   = 0;
  int maxPixel= 0;

  // starts for bounding box set, spotNo
  this(int num, int x, int y) {
    spotNo = num;
    x1 = x2 = x;           // start bounding box bounding one pixel
    y1 = y2 = y;
    sum = pixels = sumxs = sumys = maxPixel = 0;
    this.x = this.y = 0.0; 
    return;
  }

  string spotTitles(int maxval) {
     auto str = std.string.format("spotNo,x,y,sum,pixels,maxPixel(%d),x1,y1,x2,y2",maxval);
     return str;
  }

  string spotDataLine() {
      auto line = std.string.format("%d,%g,%g,%d,%d,%d,%d,%d,%d,%d",
          spotNo, x, y, sum, pixels, maxPixel, x1, y1, x2, y2);
      return line;
  }

  void addPixel(int xloc, int yloc, int value) {
     // expand bounding box if warranted
     if (x1 > xloc) { x1 = xloc; }
     if (y1 > yloc) { y1 = yloc; }
     if (x2 < xloc) { x2 = xloc; }
     if (y2 < yloc) { y2 = yloc; }
     //
     if (value > maxPixel) { maxPixel = value; }
     pixels += 1;
     sum    += value;
     sumxs  += xloc;
     sumys  += yloc;
     return;
  }

  // called when spot is done.  
  // calculates centroid, currently not weighed by intensity
  void calcVals() {     
     if (pixels > 0) {
        x = cast(double) sumxs / cast(double) pixels;
        y = cast(double) sumys / cast(double) pixels;        
     }
  }
}


/*
 * A new Spotter is created with the image to be analyzed, and data
 * on which spots to keep and delete (currently max size and min size).
 * This is based on the Fred Maxfield program kspots (Ken [Dunn] spots)
 * which was written by Mike Hillmeyer.  The program was designed to
 * trim spots so that the boundaries were full-width-half-maximum.  The
 * assumption was that they were diffraction-limited endosomes taken
 * on a conventional fluorescence scope where FWHM gives a measure of
 * fluorescence robust against small focus shifts.
 * I don't know (and never knew) the algorithm used in the original VAX
 * version to find contiguous regions above a threshold, but the method
 * here of writing the spot number to a scratch image for each pixel in each 
 * spot was present in the original program, and allows some tricks to operate
 * on pixels after the spots are isolated.  This version also keeps track of a bounding
 * rectangle surrounding each spot in case more involved image analysis needs to be
 * applied to individual spots.
 *
 * A queue-based algorithm is used to find contiguous regions.  It may not be the fastest
 * but is very simple to implement, and unlikely to blow up.  Recursion is also simple,
 * faster, but could blow the stack for very large contiguous regions.  The queue approach
 * searches an expanding boundary (breadth first) while the recursive alternative 
 * could load all pixels on the stack prior to unloading and dealing with them.  
 * This version searches four-connectivity, i.e., pixels are considered contiguous if 
 * they are directly above, below, left, right.  Diagonals are not considered, but 
 * would be easy to add if desired.  Faster, but more complex algorithms exist but
 * aren't considered worth the effort by me.  Others may disagree.
 * 
 * Right now, saturated pixels ARE NOT FLAGGED.  This should be fixed if the 
 * program is used for ratio imaging. 
 * 
 * Spotter should be initialized, and normally the method kspots() run, then
 * printSpots(filename) called if desired.  The individual subroutines called by kspots() 
 * for obtaining, trimming or deleting spots can also be used separately
 * spots() writes on scratch image, obtains data on the objects
 * trimSpots(fraction) zeroes pixels in each spot below the desired fraction of the 
 * maximum pixel.  
 * printSpots(filename) writes a report on all spots isolated
 * 
 */

class Spotter {
   enum MAXSPOTS=100_000;

   // set on initialization
   int threshold;
   int maxSize;
   int minSize;
   double fraction;
   Iplane img;

   // intermediate values
   SpotQueue q;
   Iplane scratch;

   // image characteristics
   int xs, ys, maxval;

   // output information
   Iplane output;
   SpotInfo[MAXSPOTS] spots;
   int numSpots; 

   // clear spots data, does not alter state of "output" or image size, maxval info
   // used after run of "spots" if to be rerun
   void resetData() {
       int i;
       q.reset();
       //writeln("Queue reset.");
       scratch.img[]    = 0;
       //writeln("Scratch image zeroed.");
       spots[]          = SpotInfo(-1, -1, -1);
       numSpots         = 1;   
   }  

   this(Iplane imgIn, int t, int max, int min, double fract) {
      int i;
      if (imgIn.tuple != 1) {
         writefln("imgIn.tuple = %d", imgIn.tuple);
         throw new Exception("Multi-tuple image fed to Spotter.  Explode it first!");
      }
      img        = imgIn;  
      xs         = imgIn.xsize;
      ys         = imgIn.ysize;
      maxval     = imgIn.maxval;
      threshold  = t;
      maxSize    = max;
      minSize    = min;
      fraction   = fract;
      scratch    = new Iplane("scratch",xs,ys,1,maxval); // ignore maxval  
      output     = new Iplane("spotted",xs,ys,1,maxval);
      for (i=0; i<output.img.length; i++) {
         output.img[i] = img.img[i];
      }
      numSpots   = 1;   // written to scratch image to mark searched, so zero leads to infinite loop
   }

   int threshPixel(int p) {
      return (p>=threshold)?p:0;
   }

   void spot() {
      int row, col, index, i;
      // zero below threshold
      //output   = img.applyPixel(&this.threshPixel);
      for (i=0; i<xs*ys; i++) {
         if (output.img[i] < threshold) {
            output.img[i] = 0;
         }
      }
     //writeln("Output image thresholded.");
     // loop through image to do spotting 
     for (row = 0; row < ys; row++) {
       for (col = 0; col < xs; col++) {
            index = row * xs + col;
            if ((output.img[index] != 0) && (scratch.img[index] == 0)) {
              SpotInfo spot = doSpot(row, col, numSpots);
	      spots[numSpots] = spot;
              numSpots += 1;
            }
          }
     }
   }

   SpotInfo doSpot(int row, int col, int spotNumber) {

      SpotInfo data;
      int i, j;                  /* utility variables for loops */
      Point xy;                  /* location of current spot */
      int pixVal;                /* value of pixel being queried */
      bool started = false;

      q.reset();
      xy.x = col;
      xy.y = row;

      q.enq(xy);                 /* always start with current pixel */
      scratch.img[xy.y*xs+xy.x] = spotNumber;

      /* first find the block of contiguous pixels and gather statistics */
     while ( ! q.empty) {
        xy = q.deq();
        if (!started) {
           data = SpotInfo(spotNumber, xy.x, xy.y);
           started = true;
        }
        data.addPixel(xy.x, xy.y, img.img[xy.y*xs+xy.x]);
        enqueueNeighbors(xy.y,xy.x,spotNumber);
     }
     data.calcVals();           // finish job on data
     return data;
   }

  
   void enqueueNeighbors(int row, int col, int spotNo) {
      Point xy;                  /* location of current spot */
      
      scratch.img[row*xs+col] = spotNo;      

      if ((col > 0) && (output.img[row*xs+col-1] != 0) && 
	    (scratch.img[row*xs+col-1] == 0)) { 
         scratch.img[row*xs+col-1] = spotNo;
         xy.x = col - 1;
         xy.y = row;
         q.enq(xy);
      }

      if ((col < (xs-1)) && (output.img[row*xs+col+1] != 0) && 
	   (scratch.img[row*xs+col+1] == 0)) { 
         scratch.img[row*xs+col+1] = spotNo;
         xy.x = col + 1;
         xy.y = row;
         q.enq(xy);
      }

      if ((row > 0) && (output.img[(row-1)*xs+col] != 0) && 
	    (scratch.img[(row-1)*xs+col] == 0)) { 
          scratch.img[(row-1)*xs+col] = spotNo;
          xy.x = col;
          xy.y = row - 1;
          q.enq(xy);
      }

      if ((row < (output.ysize-1)) && (output.img[(row+1)*output.xsize+col] != 0) && 
	   (scratch.img[(row+1)*xs+col] == 0)) { 
         scratch.img[(row+1)*xs+col] = spotNo;
         xy.x = col;
         xy.y = row + 1;
         q.enq(xy);
      }
      return;
   }

   // spots[numSpots]
   //
   // spots has to be run first
   void trimSpots(double fract) {
      int i, row, col;
      double max;
      for (row=0; row<ys; row++) {
         for (col=0; col<xs; col++) {
              int temp = scratch.img[row*scratch.xsize+col];
              if (temp != 0) {
                 double m = cast(double) output.img[row*output.xsize+col];
                 max = cast(double) spots[temp].maxPixel;
                 // zero pixels sufficiently below brightest pixel in spot
                 if (m < (max * fract)) {
                    output.img[row*output.xsize+col] = 0;
                    scratch.img[row*xs+col]          = 0;
                 }
              }
          } 
      }
   }


   Iplane fastSpots() {
      this.spot();
      this.trimSpots(fraction);
      this.resetData();
      this.spot();
      deleteSpots(minSize, maxSize);
      this.resetData();
      this.spot();
      return output;
   }

   Iplane kspots() {
       double tempFract;
       //writeln("Kspots entered.");
       double increment = fraction / 5.0;
       for (tempFract = increment; tempFract < fraction; tempFract+= increment) {
           //writefln("tempFract = %g", tempFract);
           this.resetData();
           //writeln("Data reset.");
           this.spot();
           //writeln("Image spotted.");
           this.trimSpots(tempFract);
           //writeln("Spots trimmed.");
       }
       //writefln("Main spot loop done. %d spots identified", numSpots );
       this.resetData();            
       this.spot();
       this.trimSpots(fraction);
       this.resetData();
       this.spot();
       //writefln("minSize=%d;  maxSize=%d  DELETING . . .", minSize, maxSize);
       deleteSpots(minSize, maxSize);  // kill below minimum or above maximum
       this.resetData();
       this.spot();                   // final spotting to produce data
       //writefln("Complete spotting done. %d spots identified", numSpots );
       return output;
   }

/* 
 * Scratch image allows matching pixels with spots in spots[i]
 * [i] is written on the scratch image for each pixel in spot[i]
 * This allows easy removal of a spot if its spot number is known.
 * In current version, number of pixels is tested and spot deleted
 * above a maximum value or below a minimum value, but any function that
 * can make a decision based on data in SpotInfo, and on the output and 
 * scratch image could be used.
 */
   void deleteSpots(int min, int max) {
      int i;
      int deletedPixels = 0;
      bool[] killList = new bool[numSpots];
      killList[] = false;
      // careful . . . numSpots[0] is never  used
      for (i=1; i<numSpots; i++) {
         int pixels = spots[i].pixels;
         if ((pixels>max)||(pixels<min)) {
            killList[i] = true;
         }         
      }
      for (i=0; i<(xs*ys); i++) {
         if (killList[scratch.img[i]]) {
            ++deletedPixels;
            scratch.img[i] = 0;
            output.img[i]  = 0;
         }
      }
      //writefln("%d deleted pixels.", deletedPixels);
   }

   // deleteSpot not as efficient as doing them in a batch
   // but easy to base on complex criteria
 
  void deleteSpot(Iplane deleted, int spotNo) {
      int row, col, index;
      if ((deleted.xsize != output.xsize) ||
          (deleted.ysize != output.ysize) ||
          (deleted.tuple != output.tuple))     {
         throw new Exception("spots.deleteSpot() mismatched deleted size");
      }
      SpotInfo s = spots[spotNo];
      for (row=s.y1; row<=s.y2; row++) {
          for (col=s.x1; col<=s.x2; col++) {
             index = row*output.xsize+col;
             if (scratch.img[index] == spotNo) {
                deleted.img[index] = output.img[index];
                output.img[index]  = 0;
                scratch.img[index] = 0;
             }
          }
      }
   }

/*

   Iplane deleteSpotsReturningDeleted(bool function killTest(Iplane,Iplane,spotInfo)) {
      Iplane deleted = output.clone();
      deleted.img[]  = 0;
      int i;
      int deletedPixels = 0;
      bool[] killList = new bool[numSpots];
      killList[] = false;
      // careful . . . numSpots[0] is never  used
      for (i=1; i<numSpots; i++) {
         int pixels = spots[i].pixels;
         if (killTest(output, scratch, spots[i]) {
            killList[i] = true;
         }         
      }
      for (i=0; i<(xs*ys); i++) {
         if (killList[scratch.img[i]]) {
            ++deletedPixels;
            deleted.img[i] = output.img[i];
            scratch.img[i] = 0;
            output.img[i]  = 0;
         }
      }
      return deleted;
   }

*/
   // *** printSpots

  // Point[] extractEdges(spotNo) {
  //    
  // }

   void printSpotTitles() {
       writef("spotNo,x,y,sum,pixels,maxPixel(%d),x1,y1,x2,y2",img.maxval);
   }

   void printSpots(int start, int finish) {
      printSpotTitles();
      foreach (s; spots[start..finish]) {
         writef("%d,%g,%g,%d,--%d--,%d,%d,%d,%d,%d",s.spotNo,s.x,s.y,s.sum,s.pixels,s.maxPixel,
                   s.x1, s.y1, s.x2, s.y2);
      }
   }

  string spotTitles() {
     auto str = std.string.format("spotNo,x,y,sum,pixels,maxPixel(%d),x1,y1,x2,y2",img.maxval);
     return str;
  }

  string spotDataLine(int spotNo) {
      auto s = spots[spotNo];
      auto line = std.string.format("%d,%g,%g,%d,%d,%d,%d,%d,%d,%d",
          s.spotNo, s.x,s.y,s.sum,s.pixels,s.maxPixel,s.x1, 
          s.y1, s.x2, s.y2);
      return line;
  }

  string allSpotData(int start, int finish) {
      //SpotInfo s;
      auto output = "";
      string st = spotTitles() ~ "\n";
      foreach (s; spots[start..finish]) {
         auto line = std.string.format("%d,%g,%g,%d,%d,%d,%d,%d,%d,%d\n",s.spotNo,
            s.x,s.y,s.sum,s.pixels,s.maxPixel,s.x1, s.y1, s.x2, s.y2);
         output ~= line;
      }
      return st ~ output;
  }

  string allSpotData() {
    return allSpotData(1,numSpots);
  }

   void printSpots() {
      this.printSpots(0, numSpots);
   }
 }

unittest {
   writeln(" . . . Testing queueing object and functions . . .");
   Point p1, p2, p3, test;
   p1.x = 2;
   p1.y = 3;
   p2.x = 4;
   p2.y = 7;
   p3.x = 10;
   p3.y = 1; 
   SpotQueue q = SpotQueue();   
   q.enq(p1);
   test = q.deq();
   assert(test == p1);
   writeln("Single enqueue, dequeue passed."); 
   q.enq(p3);  q.enq(p2); q.enq(p1);
   test = q.deq();
   assert(test == p3);
   test = q.deq();
   assert(test == p2);
   test = q.deq();
   assert(test == p1);
   writeln("Triple enqueue, dequeue passed.");
   int i;
   // force wrap-around
   for (i=0; i < q.MAXQUEUE+10; i++) {
      q.enq(p1); q.enq(p2);
      test = q.deq();
      test = q.deq();
   }
   q.enq(p3); q.enq(p2); q.enq(p1);
   test = q.deq();
   assert(test==p3);
   test = q.deq();
   assert(test == p2);
   test = q.deq();
   assert(test == p1);
   writeln("Test of wrap-around queueing passed.");
   writeln("All queueing tests passed.");
   writeln("Now testing spotter.");
   writeln("No spotter tests yet!!!!!!!!!!!!!!!!!!!");
}


