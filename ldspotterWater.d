
import std.stdio;
import std.c.stdlib;
import std.getopt;
import std.string;
import std.math;
import std.random;
import std.conv;
import image;
import spots;
import convolve;
import circlefit;
import watershedLib;
import wrapper;

enum Color { uncolored, red, yellow, green, blue }

struct ProcessResults {
   Iplane processed;
   Iplane accepted;
   Iplane deleted;
   SpotInfo[] spotList;
   CircData[] circles;
   SpathSolver[] solvers;
   Color[] color;
   int numSpots;
}
 
// Bresenham Algorithm for drawing circle
// setPixel clips but if circle is super-large a lot of time could be spent
// trying to draw off-screen

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
   auto line = header ~ "," ~ g.spotDataLine(1);  
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
   const maxIterations  = 100;
   const double epsilon = 0.1;      // 10% of pixel is good enough fit for govmut work
   double ssqTracked,oldSsqTracked;  // tracked sum of squares
   int i;
   double n;
   Iplane img = o.output;
   ColorOrKill result = ColorOrKill();
   
   oldSsqTracked = -50_000;         // Clearly invalid value
   for (i=0; i < maxIterations; i++) {
      s.nextIteration();
      if ((i % 10) == 0) {
        ssqTracked = s.sumOfSpathSquares();
        if ((abs(oldSsqTracked - ssqTracked)/cast(float)s.n) < epsilon) {
           writefln("%d iterations to fit circle.", i); 
           break;
        }
        oldSsqTracked = ssqTracked;
      }
      writeln("***500 iterations!!");
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

ProcessResults processDropletImage(Iplane img, File outfile, int minPixels, ColorOrKill function(OutlineGetter, CircData, ref SpathSolver, int) killTest) 
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
    threshold  = determineThreshold(img);
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
       if (colorkill.red >= 0)   { bresenhamCircle(red, at, bt, rt, colorkill.red); }
       if (colorkill.green >= 0) { bresenhamCircle(green, at, bt, rt, colorkill.green); }
       if (colorkill.blue >= 0)  { bresenhamCircle(blue, at, bt, rt, colorkill.blue); }
       if ((colorkill.red >= 0) && (colorkill.green >= 0)) {
          returnable.color[o] = Color.yellow;
       } else if (colorkill.red >= 0) {
          returnable.color[o] = Color.red;
       } else if (colorkill.green >= 0) {
          returnable.color[o] = Color.green;
       } else {
          returnable.color[o] = Color.blue;
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

// writeData should deal directly with the spotter since that knows the number of spots
// current version gets really ugly concatenated spot list 
// and number of spots == original + watershed spots
// this mess should be cleaned up.

void writeData(string fileName, SpotInfo[] si, Color[] c, SpathSolver[] ss, int numSpots, int maxval) {
   int i;
   string header, line;
   File f = File(fileName, "w");
   // the strings returned from spots and solvers have no consistent format!! fix!
   header = si[0].spotTitles(maxval)  ~ "circleColor" ~ ss[0].titlesAsString(); // titlesAsString has \n
   f.write(header);
   string colorAsString;
   // spot[0] is always empty
   for (i=1; i<numSpots; i++) {
      if (i < c.length) {
         colorAsString = to!string(c[i]);
      } else {
         colorAsString = "watershed";
      }
      if (i < ss.length) {
         line = si[i].spotDataLine() ~ "," ~ colorAsString ~ "," ~ ss[i].valuesAsString();  // vAS /n terminated
      } else {
         // the following give no info why circle is not solved
         // actually indicates that spot was isolated after watershed
         // for now, this has just to be known, but is not the right way to do things
         line = si[i].spotDataLine() ~ ",0.0,0.0,0.0,0.0,0,0.0,0.0\n";   
      } 
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
  string total      = "";   // accepted + watershed processed
  string datafile   = "";   // data to write to after watershed run
  string datafile_pre = ""; // data to write before (includes unprocessed clusters)
  int minPixels     = 6;    // minimum number of pixels to accept
  int threshold     = -1;   // negative value leads to automated determination of threshold
  real gauss        = 2.5;  // radius to blur image before applying watershed -- may have to be increased from default if massive image noise
  bool nowatershed  = false;  // skip watershed if set
}

Iplane chooseColor(Iplane img, int channel) {
  Iplane[] imgs;

  if (img.tuple > 1) {
    imgs = img.explode();
  }
  switch(channel) {
  case NO_COLOR: {
    if (img.tuple != 1) {
      writeln("Color image loaded.  MUST specify which channel with --color ## option");
      writeln("   --color 0       for red channel");
      writeln("   --color 1       for green channel");
      writeln("   --color 2       for blue channel");
      exit(1);
    }
    break;
  }
  case 0: {     // red channel
    if (img.tuple == 1) {
      writeln("  *** WARNING *** Greyscale image.  Red channel requested . . . using only channel . . . ");
      break;
    } else {
      img = imgs[0];
      break;
    }
  }
  case 1: {    // green channel
    if (img.tuple == 1) {
      writeln("  *** WARNING *** Greyscale image.  Green channel requested . . . using only channel . . . ");
      break;
    } else {
      img = imgs[1];
      break;
    }
  }
  case 2: {    // blue channel
    if (img.tuple == 1) {
      writeln("  *** WARNING *** Greyscale image.  Blue channel requested . . . using only channel . . . ");
      break;
    } else if (img.tuple < 3) {
      writeln(" *** Blue image (channel 3) requested, but image doesn't have three channels.");
      writeln(" *** Program terminating.");
      exit(1);
    } else {
      img = imgs[2];
      break;
    }
  }
  default: {
    writeln("Invalid channel requested.  --color 0 | --color 1 | --color 2 | --color -1 are only valid options. --color -1 for greyscale");
    writeln(" *** Program terminating");  
    exit(1);
  }
  }
  return img;
}

void main(string[] args) {
  Iplane color, accepted, killed;
  OutlineGetter outliner;
  Iplane processed;
  Point []outlineArray;
  CircData c;
  SpathSolver s;
  string input, circles, output, rejected, data;



  Parms p;
  getopt(args,
	 "in", &p.inImg,                   // input lipid droplet greyscale image
         "color", &p.color,
         "circles", &p.circlesImg,         // display image showing color-coded circles around all found objects
         "accepted", &p.outImg,                 // accepted (acceptably circular) structures
         "rejected", &p.rejectedImg,       // rejected (non-circular) structures 
         "data", &p.datafile,              // output file for data
         "predata", &p.datafile_pre,
         "minsize", &p.minPixels,          // smallest object accepted  
         "thresh", &p.threshold,           // threshold, if you don't want to automatically set
         "gauss", &p.gauss,                // default is usually okay, increase if image too noisy and makes watershed go crazy
         "total", &p.total,
         "nowatershed", &p.nowatershed,    // skip watershed
	 );
  if (p.inImg == "") { 
     writeln("Usage ldspotterWater --in inputImg.pgm [--color 2] [--minsize 6 --thresh -1 --circles circles.ppm --gauss 2.5 --accepted outputimg.pgm --rejected rejectedspots.pgm --total total.pgm --predata nowatershed.dat --data datafile.dat]");
     writeln("Input image not specified.  Program terminating.");
     exit(1);
  } 

  input    = p.inImg;
  circles  = p.circlesImg;
  output   = p.outImg;
  rejected = p.rejectedImg;
  data     = p.datafile;
  auto data_pre = p.datafile_pre;

  int minPixels = p.minPixels;
  // do test on single lipid droplets
  Iplane img    = readPNM(input);
  img           = chooseColor(img, p.color);
  File outfile  = File("temp.dat", "w");         // *** shouldn't have to do, clean up later
  auto results  = processDropletImage(img, outfile, minPixels, &killTest);   // innards totally different from ldspotter!! Redo both later
  color         = results.processed;
  accepted      = results.accepted;
  killed        = results.deleted;
  if (data_pre != "") {
      writeData(data_pre, results.spotList, results.color, results.solvers, results.numSpots, img.maxval);
  }
  if (p.nowatershed) {    // if watershed not run, finish up and quit
    if (data != "") {
      writeData(data_pre, results.spotList, results.color, results.solvers, results.numSpots, img.maxval);
    }
    if (circles != "") {
      writePNM(circles, color);
    }
    if (output != "") {
      writePNM(output, accepted);
    }
    if (rejected != "") {
      writePNM(rejected, killed);
    }
    exit(0);
  }
  auto w        = Watershedder(killed.gauss(p.gauss).invert());
  w.process();
  w.errorCount();   
  Iplane waterlines = w.waterLines(killed);
  writePNM("waterlines.pgm", waterlines);
  Iplane sum        = accepted.clone();
  sum.img[]        += waterlines.img[];                // may need to adjust maxval
  auto wSpot        = new OutlineGetter(sum, 1, 1_000_000, 6, 0.5);
  wSpot.fastSpots();

  if ((data != "") && (data !is null)) {
     //writeData(data, results.spotList, results.color, results.solvers, numSpots, img.maxval);
     string sd = wSpot.allSpotData();
     auto d = File(data,"w");
     d.writefln("%s", sd);
     d.close();
  }
  if (circles != "") {
     writePNM(circles, color);
  }
  if (output != "") {
    writePNM(output, accepted);
  }
  if (rejected != "") {
    writePNM(rejected, killed);
  }
  if (p.total != "") {
    writePNM(p.total, sum);
  }
}






 