
import std.stdio;
import std.c.stdlib;
import std.string;
import std.math;
import std.getopt;
import image;
import spots;

struct Parms {
    string img1Name;    // should have been median filtered already
    string img2Name;
    string data1Name;   // data files to write
    string data2Name;
    int minPixels  = 6;
    int maxPixels  = 1_000_000;
    int thresh1    = 15;
    int thresh2    = 15;
    double fract1  = 0.5;       // spot fraction img1 default to FWHM
    double fract2  = 0.5;       // spot fraction img2 default to FWHM
    double coloc1  = 0.3;       // fraction colocalization first must have with second
    double coloc2  = 0.3;       // fraction colocalization second must have with first
}

struct ColocResults {
   SpotInfo spot1;
   int pixelsHit;
   int pixelsMissed;
   double fractColoc = 0.0;
}

class Colocalizer {
   Spotter s1, s2;
   ColocResults[] forward_12;
   ColocResults[] backward_21;
   double coloc1, coloc2;             // fraction colocalization threshold
   bool memoized = false;
   Iplane coloc1Yes2, coloc1No2, coloc2Yes1, coloc2No1;

   this(Spotter s1, Spotter s2, double coloc1, double coloc2) {
      this.s1          = s1;
      this.s2          = s2;
      this.coloc1      = coloc1;
      this.coloc2      = coloc2;
      this.forward_12  = new ColocResults[s1.numSpots];
      this.backward_21 = new ColocResults[s2.numSpots];
   }

   private void forwardColoc(Spotter s1, Spotter s2, double colocThresh,  out Iplane kept, 
             out Iplane deleted, out ColocResults[] data) {

      //int i, pix;
      writeln("In forwardColoc");
      data          = new ColocResults[s1.numSpots+1];
      //writefln("About to do forward coloc spot copy on %d spots.", s1.numSpots);
      //writefln("Length of s1 is %d.", s1.spots.length);
      for (int i=0; i<s1.numSpots; i++) {
         data[i].spot1 = s1.spots[i];
      }
      writeln("Forward coloc spot copy done.");
      Iplane saveds = s1.scratch.clone();   // spot deletion side effects scratch image
      Iplane savedo = s1.output.clone();    // save originals so Spotter unchanged
      deleted       = s1.output.blankClone();
      writeln("About to start first foreach");
      foreach (int i, int pix; s1.scratch.img) {
         if (pix != 0) {
            if (s2.scratch.img[i] != 0) {
               ++data[pix].pixelsHit;
            } else {
               ++data[pix].pixelsMissed;
            }
         }
      }
      writeln("Finished first foreach, about to start second one.");
      foreach (i, ref ColocResults d; data) {
         double c, m;
         c = cast(double) d.pixelsHit;
         m = cast(double) d.pixelsMissed;
         d.fractColoc = c / (c + m);
         if (d.fractColoc < colocThresh) {
            s1.deleteSpot(deleted, cast(int) i);
         }
      }
      writeln("Finished second foreach");

      kept = s1.output;
      // restore old condition
      s1.output = savedo;
      s1.scratch= saveds;
      return;
   } 

   Iplane[] coloc() {
      Iplane[] answer;
      Iplane cScratch1, cScratch2;
      int numSpots1, numSpots2;
      if (memoized) {
         answer = [coloc1Yes2, coloc1No2, coloc2Yes1, coloc2No1];
         return answer;
      }
      writeln("In coloc, first calling.");
      // test that dimensions are the same
      numSpots1 = s1.numSpots;
      numSpots2 = s2.numSpots;
      cScratch1 = s1.scratch;
      cScratch2 = s2.scratch;
      writeln("About to forwardColoc");
      forwardColoc(s1, s2, coloc1, coloc1Yes2, coloc1No2, forward_12);
      writeln("Completed first forwardColoc");
      forwardColoc(s2, s1, coloc2, coloc2Yes1, coloc2No1, backward_21);
      writeln("Completed second forwardColoc");
      answer = [coloc1Yes2, coloc1No2, coloc2Yes1, coloc2No1];
      memoized = true;
      return answer;
   }

   void writeHeader(File f) {
      f.writeln("spotNum,x0,y0,pixels,sum,maxPixels,pixelsHit,pixelsMissed,fractColoc");
   }

   void writeResults(File f, bool first) {
      //int i;
      ColocResults[] r;
      r = (first)? forward_12 : backward_21;
      foreach(int i, ColocResults n; r) {
         f.writefln("%d,%g,%g,%d,%d,%d,%d,%d,%g", i, n.spot1.x, n.spot1.y, n.spot1.pixels,
            n.spot1.sum, n.spot1.maxPixel, n.pixelsHit, n.pixelsMissed, n.fractColoc);
      }
   }
}

void printUsage() {
    writeln("coloc --img1 i1.pgm --img2 i2.pgm [ --min 6 --max 1000 --thresh1 15 --thresh2 15 fract1 0.5 --fract2 0.5 --coloc1  0.3 --coloc2 0.3 --data1 dataForward.dat --data2 dataReverse.dat]");
}

void main(string[] args) {
    Parms parms;
    Iplane img1, img2, proc1, proc2, scratch1, scratch2;
    Spotter spotter1, spotter2;
    try {
       getopt(args,
          "img1",    &parms.img1Name,
          "img2",    &parms.img2Name,
          "data1",   &parms.data1Name,
          "data2",   &parms.data2Name,
          "min",     &parms.minPixels,
          "max",     &parms.maxPixels,
          "thresh1", &parms.thresh1,
          "thresh2", &parms.thresh2,
          "fract1",  &parms.fract1,
          "fract2",  &parms.fract2,
          "coloc1",  &parms.coloc1,
          "coloc2",  &parms.coloc2
       );
       if ((parms.img1Name=="")||(parms.img2Name=="")) {
          writeln("Two valid input images required but not obtained.");
          throw new Exception("Not two input images.");
       } else {
          img1 = readPNM(parms.img1Name);
          img2 = readPNM(parms.img2Name);
       }
    } catch (Exception e) {
       printUsage();
       writeln("Program exiting.");
       exit(1);
    }
    writeln("About to make spotters.");
    spotter1 = new Spotter(img1,parms.thresh1,parms.maxPixels,parms.minPixels,parms.fract1);
    writeln("Spotter1 made");
    spotter2 = new Spotter(img2,parms.thresh2,parms.maxPixels,parms.minPixels,parms.fract2);
    writeln("Spotter2 made.");
    proc1    = spotter1.kspots();
    writeln("Spotter1 kspotted.");
    proc2    = spotter2.kspots();
    writeln("Spotter2 kspotted");
    Colocalizer c = new Colocalizer(spotter1, spotter2, parms.coloc1, parms.coloc2);
    writeln("Colocalizer made.");
    Iplane[]imgs  = c.coloc();
    writeln("Coloc run.");
    // imgs are [0]c1yesc2 [1]c1noc2 [2]c2yesc1 [3]c2noc1
    // *** write them
    writePNM("coloc1_2yes.pgm", imgs[0]);
    writePNM("coloc1_2no.pgm", imgs[1]);
    writePNM("coloc2_1yes.pgm", imgs[2]);
    writePNM("coloc2_1no.pgm", imgs[3]);

    // *** write data to .dat file
    if (parms.data1Name != "") {
       File forward = File(parms.data1Name,"w");
       c.writeHeader(forward);
       c.writeResults(forward, true);
       forward.close;
    }

    if (parms.data2Name != "") {
       File reverse = File(parms.data2Name,"w");
       c.writeHeader(reverse);
       c.writeResults(reverse, false);
       reverse.close;
    }
}






 