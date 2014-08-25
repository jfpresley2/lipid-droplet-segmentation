 
import std.stdio;
import std.c.stdlib;
import std.string;
import std.math;
import std.getopt;
import image;
import spots;
import watershedLib;
import convolve;
import circleFit;
import wrapper;

struct ReturnableImages {
   Iplane accepted;
   Iplane rejected;
   Iplane combined;
   Iplane circled;
   Iplane colocalized;
   Iplane notColocalized;

   void writeImages(string rootname) {
      if (accepted !is null) {
         writePNM(rootname ~ "_accepted.pgm", accepted);
      }
      if (rejected !is null) {
         writePNM(rootname ~ "_rejected.pgm", rejected);
      }
      if (combined !is null) {
         writePNM(rootname ~ "_combined.pgm", combined);
      }
      if (circled !is null) {
         writePNM(rootname ~ "_circled.ppm", circled);
      }
      if (colocalized !is null) {
         writePNM(rootname ~ "_colocalized.pgm", colocalized);
      }
      if (notColocalized !is null) {
         writePNM(rootname  ~ "_notColocalized.pgm", notColocalized);
      }
   }
}

struct Parms {
    string img1Name;    // should have been median filtered already
    string img2Name;
    string color="";    // optional, can contain any two-letter combination of rgb
                        // !compatible with use of img2Name, both channels from img1
    string data1Name;   // data files to write
    string data2Name;
    int minPixels  = 6;
    int maxPixels  = 1_000_000;
    int thresh1    = -1;
    int thresh2    = -1;
    double fract1  = 0.5;       // spot fraction img1 default to FWHM
    double fract2  = 0.5;       // spot fraction img2 default to FWHM
    double coloc1  = 0.3;       // fraction colocalization first must have with second
    double coloc2  = 0.3;       // fraction colocalization second must have with first
    string returnRoot1 = "";    // I'm being lazy about this.  Return all possible
    string returnRoot2 = "";    // images or none.
    bool noWatershed1  = false;   // set to skip watershed on either image
    bool noWatershed2  = false;
    double gauss1  = 2.5;         // gauss for watershed
    double gauss2 = 2.5;
}

struct ColocResults {
   SpotInfo spot1;
   int pixelsHit;
   int pixelsMissed;
   double fractColoc = 0.0;
}

struct ProcessResults {
   Iplane processed;
   Iplane accepted;
   Iplane deleted;
   SpotInfo[] spotList;
   CircData[] circles;
   SpathSolver[] solvers;

}

void bresenhamCircle(Iplane img, int x0, int y0, int r, int fill) {
   if (img.tuple != 1) {
      throw new Exception("drawCircle can't handle multi-tuple images.");
   }
   int xs= img.xsize;
   int x = r;
   int y = 0;
   int offRadius = 1-x;
   while (x>=y) {
      img.setPixel(x+x0, y+y0, fill);
      img.setPixel(x0+y, x+y0, fill);
      img.setPixel(x0-x, y+y0, fill);
      img.setPixel(x0-y, x+y0, fill);
      img.setPixel(y+x0, y0-x, fill); 
      img.setPixel(x0-x, y0-y, fill);
      img.setPixel(x0-y, y0-x, fill);
      img.setPixel(x0+x, y0-y, fill);          

      y++;
      if(offRadius<0) {
        offRadius += 2*y+1;
      } else {
        x--;
        offRadius+=2*(y-x+1);
      }
   }
   return;
}

class Colocalizer {
   OutlineGetter s1, s2;
   ColocResults[] forward_12;
   ColocResults[] backward_21;
   double coloc1, coloc2;             // fraction colocalization threshold
   bool memoized = false;
   Iplane coloc1Yes2, coloc1No2, coloc2Yes1, coloc2No1;

   this(OutlineGetter s1, OutlineGetter s2, double coloc1, double coloc2) {
      this.s1          = s1;
      this.s2          = s2;
      this.coloc1      = coloc1;
      this.coloc2      = coloc2;
      this.forward_12  = new ColocResults[s1.numSpots];
      this.backward_21 = new ColocResults[s2.numSpots];
   }

   private void forwardColoc(OutlineGetter s1, OutlineGetter s2, double colocThresh,  
     out Iplane kept, out Iplane deleted, out ColocResults[] data) {

      //int i, pix;
      data          = new ColocResults[s1.numSpots+1];
      //writefln("About to do forward coloc spot copy on %d spots.", s1.numSpots);
      //writefln("Length of s1 is %d.", s1.spots.length);
      for (int i=0; i<s1.numSpots; i++) {
         data[i].spot1 = s1.spots[i];
      }
      Iplane saveds = s1.scratch.clone();   // spot deletion side effects scratch image
      Iplane savedo = s1.output.clone();    // save originals so Spotter unchanged
      deleted       = s1.output.blankClone();
      foreach (int i, int pix; s1.scratch.img) {
         if (pix != 0) {
            if (s2.scratch.img[i] != 0) {
               ++data[pix].pixelsHit;
            } else {
               ++data[pix].pixelsMissed;
            }
         }
      }
      foreach (i, ref ColocResults d; data) {
         double c, m;
         c = cast(double) d.pixelsHit;
         m = cast(double) d.pixelsMissed;
         d.fractColoc = c / (c + m);
         if (d.fractColoc < colocThresh) {
            s1.deleteSpot(deleted, cast(int) i);
         }
      }

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
      // test that dimensions are the same
      numSpots1 = s1.numSpots;
      numSpots2 = s2.numSpots;
      cScratch1 = s1.scratch;
      cScratch2 = s2.scratch;
      forwardColoc(s1, s2, coloc1, coloc1Yes2, coloc1No2, forward_12);
      writeln("Completed first forwardColoc");
      forwardColoc(s2, s1, coloc2, coloc2Yes1, coloc2No1, backward_21);
      writeln("Completed second forwardColoc");
      answer = [coloc1Yes2, coloc1No2, coloc2Yes1, coloc2No1];
      memoized = true;
      return answer;
   }

   void writeHeader(File f) {
      f.writeln("spotNum,x0,y0,pixels,sum,maxPixels,pixelsHit,pixelsMissed,fractColoc,tag");
   }

   void writeResults(File f, string tag, bool first) {
      //int i;
      ColocResults[] r;
      r = (first)? forward_12 : backward_21;
      // r[0] not used, r[$] not used
      foreach(int i, ColocResults n; r[1..($-1)]) {
         f.writefln("%d,%g,%g,%d,%d,%d,%d,%d,%g, %s", i, n.spot1.x, n.spot1.y, n.spot1.pixels,
            n.spot1.sum, n.spot1.maxPixel, n.pixelsHit, n.pixelsMissed, n.fractColoc, tag);
      }
   }
}

struct ColorOrKill {
   bool killed = false;
   int red     = -1;
   int green   = -1;
   int blue    = -1;
}

double diff(double a, double b) {
   double d = a - b;
   return (d >= 0.0)? d : -d;
}

ColorOrKill killTest(OutlineGetter o, CircData c, ref SpathSolver s, int spotNum) {
   int i;
   double n;
   Iplane img = o.output;
   ColorOrKill result = ColorOrKill();
   for (i=0; i<500; i++) {
      s.nextIteration();
   } 
   double ss       = s.sumOfSquares();
   n               = cast(double) s.n;
   double miss     = sqrt(ss/n) / s.rt;
   double areaCirc = 3.1415926 * s.rt * s.rt;
   double areaPix  = cast(double) o.spots[spotNum].pixels;
   double ratio    = areaCirc / areaPix;

   // use area of circles vs pixels to find fits that are grossly off i.e. failure

   if ((ratio > 1.5) || (ratio < 0.667)) { 
      goto failblue; 
   }

   // sd of fit worse than area spot spread over -- probably arc fit in bizarre way

   if (cast(ulong)miss > (diff(o.spots[spotNum].x1, o.spots[spotNum].x2)+
       diff(o.spots[spotNum].y1, o.spots[spotNum].y2))) {
    failblue:
      //writeln("Failed blue.");
      result.blue   = img.maxval; // mark fit as failed with blue (red used for usual)
      result.killed = true;
      s.rt = (s.rt > 10)? 10.0 : s.rt;  // make sure blue is visible
      s.at        = o.spots[spotNum].x; // use centroid if fit failed
      s.bt        = o.spots[spotNum].y;
      return result;
   }
   if ((miss < 0.15) && (s.rt > 5.0)) {
      result.killed = false;
      result.green  = img.maxval;
      return result;
   } else if ((( miss < 0.25) && (s.rt > 4.0) && (s.rt <= 5.0)) ||
              ((miss < 0.30) &&  (s.rt > 3.0) && (s.rt <= 4.0)) ||
              ((miss < 0.35) &&  (s.rt <= 3.0))) {
      result.killed = false;
      result.green  = img.maxval;
      result.red    = img.maxval;   // color circle yellow
      return result;
   } else {
      result.killed = true;
      result.red    = img.maxval;
      return result;
   }
   return result;
}

int determineThreshold(Iplane img) {
   return (img.maxval / 10);
}

ProcessResults processDropletImage(Iplane img1, int minPixels, int threshold, ColorOrKill function(OutlineGetter, CircData, ref SpathSolver, int) killTest)  {

    ColorOrKill colorkill;
    Iplane noisedImage, processed, accepted, deleted, writeto;
    Iplane red, green, blue;
    Iplane[3] temp;
    ProcessResults returnable;
    OutlineGetter spotter;
    int n, i, o;
    bool titleWritten = false;

    deleted       = img1.clone();
    deleted.img[] = 0;           // deleted should have same size but empty

    if (img1.tuple != 1) {
       throw new Exception("processDropletImage: can't process multi-channel droplet image.");
    }
    spotter    = new OutlineGetter(img1, threshold, img1.xsize*img1.ysize, minPixels, 0.5);
    processed  = spotter.kspots();
    red      = processed.clone();
    green    = processed.clone();
    blue     = processed.clone();
    accepted = processed.clone();

    auto spotList   = spotter.spots;
    CircData[] c    = new CircData[spotter.numSpots];
    SpathSolver[] s = new SpathSolver[spotter.numSpots];
    writefln("%d spots isolated.", spotter.numSpots);

    // *** num iterations SpathSolver fixed at 500, much less should do job
    // *** if routinely following reduction in sum of squares
    for (o=1; o<spotter.numSpots; o++) {
       c[o] = CircData(spotter.calculateOutlineArray(o));
       s[o] = SpathSolver(c[o]);
       colorkill = killTest(spotter, c[o], s[o], o);
       if (o%50==0) { writeln("On spot %d",o); }
       if (colorkill.killed) {
         spotter.deleteSpot(deleted,o);
       } else {
       }
       int at = cast(int) s[o].at;
       int bt = cast(int) s[o].bt;
       int rt = cast(int) floor(s[o].rt + 0.5);
       if (colorkill.red >= 0)   { bresenhamCircle(red, at, bt, rt, colorkill.red); }
       if (colorkill.green >= 0) { bresenhamCircle(green, at, bt, rt, colorkill.green); }
       if (colorkill.blue >= 0)  { bresenhamCircle(blue, at, bt, rt, colorkill.blue); }
    }
    writeln("Writing circles.");
    temp = [red, green, blue];
    processed = compact(temp);
    writeln("Compacted.");
    accepted = spotter.output;
    returnable.processed = processed;
    returnable.accepted  = accepted;
    returnable.deleted    = deleted;
    returnable.spotList   = spotList;
    returnable.circles    = c;
    returnable.solvers    = s;
    return returnable;
}

bool validColor(string ch, int chans) {
   if ((ch!="r")&&(ch!="g")&&(ch!="b")) {
      return false;
   }
   if ((ch=="r")&&(chans<1)) {
      return false;
   }
   if ((ch=="g")&&(chans<2)) {
      return false;
   }
   if ((ch=="b")&&(chans<3)) {
      return false;
   }
   return true;
}

void interpretColor(ref Iplane img1, out Iplane img2, string color) {
   Iplane[] imgs;
   string i1, i2;
   if (img1.tuple < 2) {
     throw new Exception("Greyscale image loaded, but color channels demanded.");
   }
   if (color.length < 2) {
     throw new Exception("Color string has < 2 chars.  Not enough information");
   }
   imgs = img1.explode();
   i1 = color[0..1];
   i2 = color[1..2];
   if (validColor(i1,img1.tuple) && validColor(i2,img1.tuple)) {
      switch(i1) {
         case "r":
            img1 = imgs[0];
            break;
         case "g":
            img1 = imgs[1];
            break;
         case "b": 
            img1 = imgs[2];
            break;
         default:
            throw new Exception("Invalid color specifying char");
      }
      switch(i2) {
         case "r":
            img2 = imgs[0];
            break;
         case "g":
            img2 = imgs[1];
            break;
         case "b": 
            img2 = imgs[2];
            break;
         default:
            throw new Exception("Invalid color specifying char");
      }
   } else {
     throw new Exception("color string must have only 'r', 'g', or 'b' in first two positions");
   }  

}

void printUsage() {
    writeln("colocWater --img1 i1.pgm --img2 i2.pgm    // mandatory parameters");
    writeln("    --color rg      // can  use instead of img2, rgb are possibles from img1");
    writeln("    // optional parameters with defaults follow");
    writeln("    --min 6 --max 1000 ");
    writeln("    --thresh1 15 --thresh2 15 fract1 0.5 --fract2 0.5");
    writeln("    --coloc1 0.3 --coloc2 0.3");
    writeln("    // optional output files, nothing written if not specified"); 
    writeln("    --data1 dataForward.dat --data2 dataReverse.dat");
    writeln("    // optional root output image filenames, a whole glob of each if specified");  
    writeln("    --returnRoot1 root --returnRoot2 root2");
}

void main(string[] args) {
    Parms parms;
    ReturnableImages output1, output2;
    Iplane img1, img2, proc1, proc2, scratch1, scratch2;
    OutlineGetter spotter1, spotter2;
    try {
       getopt(args,
          "img1",    &parms.img1Name,
          "img2",    &parms.img2Name,
          "color",   &parms.color,
          "data1",   &parms.data1Name,
          "data2",   &parms.data2Name,
          "min",     &parms.minPixels,
          "max",     &parms.maxPixels,
          "thresh1", &parms.thresh1,
          "thresh2", &parms.thresh2,
          "fract1",  &parms.fract1,
          "fract2",  &parms.fract2,
          "coloc1",  &parms.coloc1,
          "coloc2",  &parms.coloc2,
          "returnRoot1", &parms.returnRoot1,
	  "returnRoot2", &parms.returnRoot2,
	  "noWatershed1", &parms.noWatershed1,
	  "noWatershed2", &parms.noWatershed2,
	  "gauss1", &parms.gauss1,
	  "gauss2", &parms.gauss2,
       );
       if ((parms.img1Name=="")||((parms.img2Name=="")&&parms.color=="")) {
          writeln("Two valid input images required but not obtained.");
          throw new Exception("Not two input images.");
       } else {
          img1 = readPNM(parms.img1Name);
          if ((parms.color=="")||(parms.color is null)) {
             img2 = readPNM(parms.img2Name);
          } else {
            writeln("Interpreting color.");
            interpretColor(img1, img2, parms.color);
            writeln("Color interpreted.");
          }
       }
    } catch (Exception e) {
       printUsage();
       writeln("Program exiting.");
       exit(1);
    }
    if (parms.thresh1 == -1) {
      parms.thresh1 = determineThreshold(img1);
    }
    if (parms.thresh2 == 01) {
      parms.thresh2 = determineThreshold(img2);
    }

    // processDropletImage gives you a processed[0], accepted[1] and a rejected[2] image
    // run watershed on the rejected image after gaussing 2.5
    auto img1Results = processDropletImage(img1, parms.minPixels, parms.thresh1, &killTest);
    auto img2Results = processDropletImage(img2, parms.minPixels, parms.thresh2, &killTest);
    writeln("Initial processing done.");
    output1.accepted = img1Results.accepted;
    output1.circled  = img1Results.processed;
    output2.accepted = img2Results.accepted;
    output2.circled  = img2Results.processed;    
    
    Iplane i1rej_temp;
    Iplane i2rej_temp;

    if (!parms.noWatershed1) {
      i1rej_temp   = gauss(img1Results.deleted, parms.gauss1);
      auto w1      = Watershedder(i1rej_temp.invert());
      w1.process();
      w1.errorCount();
      i1rej_temp   = w1.waterLines(i1rej_temp);
    }  else {
      i1rej_temp = img1Results.deleted;           // don't blur if not needed for watershed
    }
    if (!parms.noWatershed2) {
      i2rej_temp   = gauss(img2Results.deleted, parms.gauss2);
      auto w2      = Watershedder(i2rej_temp.invert());
      w2.process();
      w2.errorCount();
      i2rej_temp   = w2.waterLines(i2rej_temp);
    } else {
      i2rej_temp = img2Results.deleted;          // don't blur if not needed for watershed
    }
    output1.rejected = i1rej_temp;
    output2.rejected = i2rej_temp;
    // produce for each a sum image (accepted + watershedded)
    auto sum1   = i1rej_temp.clone();
    auto sum2   = i2rej_temp.clone();
    sum1.img[] += img1Results.accepted.img[]; 
    sum2.img[] += img2Results.accepted.img[];
    output1.combined = sum1;
    output2.combined = sum2;

    // coloc accepted1 with sum2
    auto i1    = new OutlineGetter(img1Results.accepted, 1,1_000_000,1,0.5);
    auto i1w   = new OutlineGetter(i1rej_temp, 1,1_000_000,1,0.5);
    auto i2sum = new OutlineGetter(sum2, 1,1_000_000,1,0.5);
    i1.fastSpots();
    i1w.fastSpots();
    i2sum.fastSpots();
    writeln("Initial object isolation done.");
    // do colocalizing     
    Colocalizer f1 = new Colocalizer(i1, i2sum, parms.coloc1, parms.coloc2);
    Colocalizer fw = new Colocalizer(i1w, i2sum, parms.coloc1, parms.coloc2);
    Iplane[]f1f = f1.coloc();
    Iplane[]f1w = fw.coloc();
    output1.colocalized           = f1f[0].clone(); 
    output1.colocalized.img[]    += f1w[0].img[];
    output1.notColocalized        = f1f[1].clone;
    output1.notColocalized.img[] += f1w[1].img[];

    writeln("Forward colocalization done.");
    if (parms.data1Name != "") {
       File forward = File(parms.data1Name,"w");
       f1.writeHeader(forward);
       f1.writeResults(forward, "accepted", true);
       fw.writeResults(forward, "watershed", true);
       forward.close;
    }

    // coloc accepted2 with sum1
    writeln("Forward colocalization written.");
    auto i2    = new OutlineGetter(img2Results.accepted, 1,1_000_000,1,0.5);
    auto i2w   = new OutlineGetter(i2rej_temp, 1,1_000_000,1,0.5);
    auto i1sum = new OutlineGetter(sum1, 1,1_000_000,1,0.5);
    i2.fastSpots();
    i2w.fastSpots();
    i1sum.fastSpots();
    Colocalizer r1 = new Colocalizer(i2, i1sum, parms.coloc2, parms.coloc1);
    Colocalizer rw = new Colocalizer(i2w, i1sum, parms.coloc2, parms.coloc1);
    Iplane[]r1f    = r1.coloc();
    Iplane[]r1w    = rw.coloc();
    output2.colocalized           = r1f[0].clone(); 
    output2.colocalized.img[]    += r1w[0].img[];
    output2.notColocalized        = r1f[1].clone;
    output2.notColocalized.img[] += r1w[1].img[];
    writeln("Reverse colocalization done.");

    if (parms.data1Name != "") {
       File reverse = File(parms.data2Name,"w");
       r1.writeHeader(reverse);
       r1.writeResults(reverse, "accepted", true);
       rw.writeResults(reverse, "watershed", true);
       reverse.close;
    }
    writeln("Reverse colocalization written.");
    if (parms.returnRoot1 != "") {
       output1.writeImages(parms.returnRoot1);
    }
    if (parms.returnRoot2 != "") {
       output2.writeImages(parms.returnRoot2);
    }

    // imgs are [0]c1yesc2 [1]c1noc2 [2]c2yesc1 [3]c2noc1

}

