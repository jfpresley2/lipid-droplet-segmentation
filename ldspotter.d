
import std.stdio;
import std.c.stdlib;
import std.getopt;
import std.string;
import std.math;
import std.random;
import image;
import spots;
import convolve;
import circlefit;
import wrapper;

enum Color { uncolored, red, yellow, green, blue };

struct ProcessResults {
   Iplane processed;
   Iplane accepted;
   Iplane deleted;
   SpotInfo[] spotList;
   CircData[] circles;
   Color[] color;
   SpathSolver[] solvers;
   int numSpots;
}
 
// Bresenham Algorithm for drawing circle
// setPixel clips, so pixels outside image boundary will be ignored
// but if circle is super-large a lot of time could be spent not drawing these pixel.

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

struct CircStats {
   int sumPixels;
   int pixels;
}

void titleLine(File f, string header, OutlineGetter g, SpathSolver s) {
   f.write(header ~ g.spotTitles() ~ s.titlesAsString());
   return;
}

void printLine(File f, string header, OutlineGetter g, SpathSolver s) {
   if (g.numSpots != 2) {
      f.writefln("%s,Multiple spots",header);
      return;
   }
   auto line = header ~ g.spotDataLine(1);  
   line ~= "," ~ s.valuesAsString();
   f.write(line); 
   return;
}

//   Iplane deleteSpotsReturningDeleted(bool function killTest(Iplane,Iplane,spotInfo))

/* fit radius of 5 or less, roundness is meaningless
 * miss of 0.15 or less means good
 *
*/

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
   const maxIterations  = 500;
   const double epsilon = 0.01;      // 1% of pixel is good enough fit for govmut work
   double ssqTracked,oldSsqTracked;  // tracked sum of squares
   int i;
   double n;
   Iplane img = o.output;
   ColorOrKill result = ColorOrKill();
  
   oldSsqTracked = -50_0000;
   for (i=0; i < maxIterations; i++) {
      s.nextIteration();
      if ((i % 10) == 0) {
        ssqTracked = s.sumOfSpathSquares();
        if (abs(oldSsqTracked - ssqTracked) < epsilon) 
           break;
        oldSsqTracked = ssqTracked;
      }
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

ProcessResults processDropletImage(Iplane img, int minPixels, int thresh, ColorOrKill function(OutlineGetter, CircData, ref SpathSolver, int) killTest) 
{
    ColorOrKill colorkill;
    Iplane noisedImage, processed, accepted, deleted, writeto;
    Iplane red, green, blue;
    Iplane[] temp;
    ProcessResults returnable;
    OutlineGetter spotter;
    int n, i, o, threshold;
    bool titleWritten = false;

    deleted       = img.clone();
    deleted.img[] = 0;           // deleted should have same size but empty

    if (img.tuple != 1) {
       throw new Exception("processDropletImage: can't process multi-channel droplet image.");
    }
    if (thresh < 0) {
       threshold  = determineThreshold(img);
    } else {
      threshold  = thresh;
    }
    spotter    = new OutlineGetter(img, threshold, img.xsize*img.ysize, minPixels, 0.5);
    processed  = spotter.kspots();
    red      = processed.clone();
    green    = processed.clone();
    blue     = processed.clone();
    accepted = processed.clone();
    temp = [red, green, blue];

    auto spotList   = spotter.spots;
    CircData c[]    = new CircData[spotList.length];
    SpathSolver s[] = new SpathSolver[spotList.length];
    returnable.color = new Color[spotList.length];
    writefln("%d spots isolated.", spotter.numSpots);

    // *** num iterations SpathSolver fixed at 500, much less should do job
    // *** if routinely following reduction in sum of squares
    for (o=1; o<spotter.numSpots; o++) {
       c[o] = CircData(spotter.calculateOutlineArray(o));
       s[o] = SpathSolver(c[o]);
       colorkill = killTest(spotter, c[o], s[o], o);
       if (o%50==0) { writeln("On spot %d",o); }
       if (colorkill.killed) {
         //writeln("Spot killed");
         spotter.deleteSpot(deleted,o);
       } else {
         //writeln("Spot accepted.");
       }
       int at = cast(int) s[o].at;
       int bt = cast(int) s[o].bt;
       int rt = cast(int) floor(s[o].rt + 0.5);
       if (colorkill.red >= 0)   { 
         bresenhamCircle(red, at, bt, rt, colorkill.red);
         returnable.color[o] = Color.red;
       }
       if (colorkill.green >= 0) { 
         bresenhamCircle(green, at, bt, rt, colorkill.green);
         returnable.color[o] = Color.green; 
       }
       if (colorkill.blue >= 0)  { 
         bresenhamCircle(blue, at, bt, rt, colorkill.blue); 
         returnable.color[o] = Color.blue;
       }
       if ((colorkill.red >= 0) && (colorkill.green >= 0)) {      // overwrite for this case
	 returnable.color[o] = Color.yellow;
       } 
    }
    writeln("Writing circles.");
    processed = compact(temp);
    accepted = spotter.output;
    returnable.processed = processed;
    returnable.accepted  = accepted;
    returnable.deleted   = deleted;
    returnable.spotList  = spotList;
    returnable.circles   = c;
    returnable.solvers   = s;
    returnable.numSpots  = spotter.numSpots;
    return returnable;
}


// The combination of processDropletImage and writeData is handling things in a really ugly way
// This part of the program needs to be redone from scratch.

void writeData(string fileName, SpotInfo[] si, Color[] c, SpathSolver[] ss, int numSpots, int maxval) {
   int i;
   string header, line;
   writefln("Writing file %s.", fileName);
   File f = File(fileName, "w");
   if (si.length != ss.length) {
      writeln("spotinfo length != spathsolver length:  PROGRAMMING ERROR");
      throw new Exception("Mismatched length");
   }
   // the strings returned from spots and solvers have no consistent format!! fix!
   header = si[0].spotTitles(maxval)  ~ ",Color" ~ss[0].titlesAsString(); // titlesAsString has \n
   f.write(header);
   // spot[0] is always empty
   for (i=1; i<numSpots; i++) {
     string color;
     switch (c[i]) {
     case Color.red:    { color = "red";  break; }
     case Color.green:  { color = "green"; break; }
     case Color.blue:   { color = "blue"; break; }
     case Color.yellow: { color = "yellow"; break; }
     default:           { color = "uncolored"; break; }      // shouldn't happen
     }
      line = si[i].spotDataLine() ~ "," ~ color ~ "," ~ ss[i].valuesAsString();  // vAS /n terminated
      f.write(line);
   }
   f.close;
}

enum NO_COLOR = -1;

struct Parms {
  string inImg;             // input image
  int color = NO_COLOR;     // if not implemented, assume greyscale image
  string circlesImg = "";
  string outImg     = "";
  string rejectedImg= "";
  string datafile   = "";
  int minPixels     = 6;    // minimum number of pixels to accept
  int threshold     = -1;   // negative value leads to automated determination of threshold
}

void main(string[] args) {
  Iplane color, accepted, killed;
  Iplane[] imgs;
  OutlineGetter outliner;
  Iplane processed;
  Point []outlineArray;
  CircData c;
  SpathSolver s;
  string input, circles, output, rejected, data;
  Iplane img, colorImg;

  Parms p;
  getopt(args,
	 "in", &p.inImg,                   // input lipid droplet greyscale image
         "color", &p.color,
         "circles", &p.circlesImg,         // display image showing color-coded circles around all found objects
         "out", &p.outImg,                 // accepted (acceptably circular) structures
         "rejected", &p.rejectedImg,       // rejected (non-circular) structures 
         "data", &p.datafile,              // output file for data
         "minsize", &p.minPixels,          // smallest object accepted  
         "thresh", &p.threshold,           // threshold, if you don't want to automatically set         
	 );
  if (p.inImg == "") { 
     writeln("Usage ldspotter --in inputImg.pgm [--color 2] [--minsize 6 --threshold -1 --circles circles.ppm --out outputimg.pgm --rejected rejectedspots.pgm --data datafile.dat]");
     writeln("Input image not specified.  Program terminating.");
     exit(1);
  } 
  input    = p.inImg;
  circles  = p.circlesImg;
  output   = p.outImg;
  rejected = p.rejectedImg;
  data     = p.datafile;
  int minPixels = p.minPixels;
  // do test on single lipid droplets
  
  img    = readPNM(input);
  if (img.tuple > 1) {
    imgs = img.explode();
  }
  switch(p.color) {
  case NO_COLOR: 
    {
       if (img.tuple != 1) {
         writeln("Color image loaded.  MUST specify which channel with --color # option");
         writeln("   --color 0       for red channel");
         writeln("   --color 1       for green channel");
         writeln("   --color 2       for blue channel");
         exit(1);
       } 
       break;
    }
  case 0:     // red channel 
    {
      if (img.tuple == 1) {
        writeln("*** WARNING *** Greyscale image. Red channel requested . . . using only channel");
	break;               
      } else {
        img = imgs[0];
        break;
      }              
    }
  case 1:    // green channel
    if (img.tuple == 1) {
      writeln("*** WARNING *** Greyscale image.  Green channel requested . . . using only channel");
      break;
    } else {
      img = imgs[1];
      break;
    }
  case 2:    // blue channel
    if (img.tuple == 1) {
      writeln("*** WARNING *** Greyscale image.  Blue channel requested . . . using only channel");
      break;
    } else if (img.tuple < 3) {
      writeln(" *** Blue image (channel 3) requested, but image doesn't have 3 channels.");
      writeln(" *** Program terminating");
      exit(1);
      break;
    } else {
      img = imgs[2];
      break;
    }
  default:
    {
      writeln("Invalid channel requested.  --color 0, --color 1, --color 2 valid options");
      writeln(" *** Program terminating");
      break;
    }
  }

  string dat;
  auto results  = processDropletImage(img, minPixels, p.threshold, &killTest);
  color         = results.processed;
  accepted      = results.accepted;
  killed        = results.deleted;

  if (p.datafile != "") {
    writeln("Writing data");
    writeData(data, results.spotList, results.color, results.solvers, results.numSpots, img.maxval);
  }
  if (p.circlesImg != "") {
     writePNM(circles, color);
  }
  if (p.outImg != "") {
     writePNM(output, accepted);
  }
  if (p.rejectedImg != "") {
    writePNM(rejected, killed);
  }
}


 