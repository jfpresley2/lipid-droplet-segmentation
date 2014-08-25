import std.stdio;
import std.math;
import std.random;
import std.c.stdlib;
import std.string;
import wrapper;
import image;
import spots;


/*****

 The purpose of circle.d is to test fitting of a boundary to a circle 

 uses method of Späth (1996) Computing 57:179-185 to fit circle to data
 Claim of Späth is it can fit circle from points along an arc
 as well as spread over a complete circle

 Tests seem to indicate it works okay, even if noisy points are localized
 only on a small arc of the circle. A publication indicated it converges
 slowly vs. Levenquan-Marquandt, but I originally wanted to have code
 I completely understood in case I needed to port to GPU.  As it turns out,
 running time is just fine and doing all the fits in parallel is not needed.

*****/
  
//extern (C) int poisson(double mean);

enum double pi  = 3.1415926;
enum int Points = 15; 

struct CircData {
  double[] x;
  double[] y;

  this(double[] x, double[] y) {
     this.x = x;
     this.y = y;
  }

  this(Point[] p) {
    int i;
    x = new double[p.length];
    y = new double[p.length];
    for (i=0; i<p.length; i++) {
       x[i] = cast(double) p[i].x;
       y[i] = cast(double) p[i].y;
    }
  }

  void printVals() {
     int i;
     for (i=0; i<x.length; i++) {
        writefln("  %g,%g", x[i], y[i]);
     }
  }
}

struct Circle {
  double a;
  double b;
  double x;
  double r;
}

// fits points to a circle using the Spath algorithm
// disadvantage -- slow compared to other methods
//              -- may fail on short noisy arcs (which won't occur)
// advantages -- no library functions needed
//            -- considered very reliable
//            -- looks gpu-able
// Note -- if multiple fitting algorithms were to be used some reorganizing
//         with a common interface is worth doing.

struct SpathSolver {
  double[] xs;         // unknown points to fit
  double[] ys;
  double[] zs;         // angles, used as intermediates in fit
  int n;               // number of points to fit
  // calculated upon initializing
  double sumxs;        
  double sumys;
  // intermediate values used in calculations
  double at, bt, ct, zt, rt;
  double sumSquares = -1;
  size_t pixels;

  void printVals() {
     writefln("a=%g; b=%g, r=%g", at, bt, rt);
  }

// leading comma for title because no trailing comma for spots title
// no leading comma for values
// *** come up with better solution

  string titlesAsString() {
     return ",x0,y0,radius,sumSquares,numPixels,ss/pixels,s/radius\n";
  }

  string valuesAsString() {
     pixels = xs.length;
     if (sumSquares == -1) {
        sumSquares = this.sumOfSquares();
     }
     return std.string.format("%g,%g,%g,%g,%d,%g,%g\n",at,bt,rt,sumSquares,
          pixels, sumSquares/pixels, sqrt(sumSquares/pixels)/rt);
  }

  Circle solvedCircle() {
     Circle c = Circle();
     c.a = at;
     c.b = bt;
     c.r = rt;
     return c;
  }
  // ---- load points to fit, and do initialization ----

  this(CircData c) {
     int i;
     double e;
     this.xs = c.x;
     this.ys = c.y;
     double sumxsq, sumxy, sumysq;
     double b1xxsys, b2yxsys, b3xsys;
     sumxsq=sumxy=sumxs=sumys=sumysq=0.0;
     b1xxsys=b2yxsys=b3xsys=0.0;
     this.n = cast(int) this.xs.length;
     if (this.n != this.ys.length) {
       throw new Exception("SpathSolver given mis-matched x,y circle data");
     }
     for(i=0; i<this.n; i++) {
        sumxsq += xs[i] * xs[i];
        sumxy  += xs[i] * ys[i];
        sumxs  += xs[i];
        sumys  += ys[i];
        sumysq += ys[i] * ys[i];
        b1xxsys += xs[i] * (xs[i]*xs[i] + ys[i]*ys[i]);
        b2yxsys += ys[i] * (xs[i]*xs[i] + ys[i]*ys[i]);
        b3xsys  += xs[i] * xs[i] + ys[i] * ys[i];
     }
     // fill in matrix for solution of equation a*x = b
     // where x = [a, b, c]
     // a, b, c will be needed to determine initial guesses
     // fortran-style a(row,column) for library call
     double[] aa = new double[9];       // 3x3 matrix, access as row*3+col
     // fill first row
     aa[0*3+0] = 2.0 * sumxsq;     // 0*3.0+0
     aa[0*3+1]   = 2.0 * sumxy;
     aa[0*3+2]   = sumxs;
     // fill second row       
     aa[1*3+0]   = 2.0 * sumxy;
     aa[1*3+1]   = 2.0 * sumysq;
     aa[1*3+2]   = sumys;
     // fill third row
     aa[2*3+0]   = 2.0 * sumxs;
     aa[2*3+1]   = 2.0 * sumys;
     aa[2*3+2]   = n;
     // fill in b vector
     double[] bb = new double[3];
     bb[0]       = b1xxsys;
     bb[1]       = b2yxsys;
     bb[2]       = b3xsys;
     // put together starting values (time 0)
     // these will be updated in each iteration
     double[] xx = new double[3];
     // *** can't figure how to return an array from C without segfaulting     
     // *** gets around it for this purpose, but works only for 3x3 matrix
     solveMatrix(&xx[0], &aa[0], 3, 3, &bb[0]);
     at  = xx[0]; 
     bt  = xx[1];
     ct  = xx[2];
     rt  = sqrt(ct + at*at + bt*bt);
     zs = new double[n];
     //zs[]        = 0.0;
     return;     
  }

  void nextIteration() {
    int i;
    // sums calculated during iteration
    double sumxcosysin, sumcosz, sumsinz;
    double cosz, sinz;
    // determinant of matrix
    double d;         
    sumxcosysin = sumcosz = sumsinz = 0.0;
    // determine new rt  i.e., r(t) for next iteration
    for (i=0; i<n; i++) {
      if ((xs[i] - at) == 0.0) {
         zs[i] = pi / 2.0;
      } else {
         double temp = (ys[i] -bt) / (xs[i] - at);
         if (temp < 0.0) { temp = -temp; }
         zs[i] = atan(temp);
      }
    }
    for (i=0; i<n; i++) {
       cosz = cos(zs[i]);
       sinz = sin(zs[i]);
       if ((xs[i] - at) <= 0.0) {
          cosz = -cosz;
       }       
       if ((ys[i] - bt) <= 0.0) {
          sinz = -sinz;
       }
       sumcosz     += cosz;
       sumsinz     += sinz;
       sumxcosysin += xs[i] * cosz + ys[i] * sinz;
    }
    // calculate determinant, which must be above zero
    d = n*n * (n - (1.0/n) * ((sumcosz*sumcosz + sumsinz*sumsinz)));
    // d <= 0.0 should be very rare, just a precaution
    if (d <= 0.0) {
       throw new Exception("Below-zero determinant.  Fit failed.");
    }
    // calculate refined (new) value of rt
    rt = (n*n*sumxcosysin - n*(sumxs*sumcosz + sumys*sumsinz))/d;
    // calculate refined values for a and b
    at = (1.0/n) * (sumxs - rt * sumcosz);
    bt = (1.0/n) * (sumys - rt * sumsinz);
    return;
  }

  // distance of each point from center - radius
  double sumOfSquares() {
     int i;
     double sum = 0.0;
     for (i=0; i<n; i++) {
       sum += (rt - sqrt((at-xs[i])*(at-xs[i])+(bt-ys[i])*(bt-ys[i]))) *
              (rt - sqrt((at-xs[i])*(at-xs[i])+(bt-ys[i])*(bt-ys[i])));
     }
     sumSquares = sum;
     return sum;
  }

  // quantity minimized
  double sumOfSpathSquares() {
     int i;
     double sum = 0.0;
     for (i=0; i<n; i++) {
       sum += (ys[i] - bt - rt*sin(zs[i])*sin(zs[i]))*(ys[i] - bt - rt*sin(zs[i])*sin(zs[i]))
               + 
              (xs[i] - at - rt*cos(zs[i]))*(xs[i] - at - rt*cos(zs[i]));
     }
     sumSquares = sum;
     return sum;
  }
}

CircData makeData(double x0, double y0, double r, int points) {
   int i;
   double angle = 0.0;
   CircData c;
   c.x       = new double[points];
   c.y       = new double[points];
   auto increment = 2*pi / cast(double)(points+1);
   for (i=0; i<points; i++) {
      angle += increment;
      c.x[i] = r * cos(angle) + x0;
      c.y[i] = r * sin(angle) + y0;
   }
   return c;
}


CircData addNoise(CircData c, double sigma) {
   int i;
   int n = cast(int) c.x.length;
   double randNum;
   CircData o;
   o.x = new double[n];
   o.y = new double[n];
   for (i=0; i<n; i++) {
      randNum = gaussian(0.0, sigma);
      o.x[i] = c.x[i] + randNum;
      randNum = gaussian(0.0, sigma);
      o.y[i] = c.y[i] + randNum;
   }
   return o;
}


class OutlineGetter : Spotter {

    this(Iplane imgIn, int t, int max, int min, double fract) {
        super(imgIn,t,max,min,fract);
    }

    private bool leftIs(int row, int col, int pixVal) {
      if ((row > 0) && (scratch.img[(row-1)*xs+col] == pixVal)) { 
         return true;
      }
      return false;
    }

    private bool rightIs(int row, int col, int pixVal) {
      if ((row < (ys-1)) && (scratch.img[(row+1)*xs+col] == pixVal)) { 
         return true;
      }
      return false;
    }

    private bool upIs(int row, int col, int pixVal) {
      if ((col > 0) && (scratch.img[row*xs+col-1] == pixVal)) { 
         return true;
      }
      return false;
    }

    private bool downIs(int row, int col, int pixVal) {
      if ((col < (xs-1)) && (scratch.img[row*xs+col+1] == pixVal)) { 
         return true;
      }
      return false;
    }

    private bool isEdge(int row, int col, int spotNo) {
       if ((scratch.img[row*xs+col] == 0) && 
             (upIs(row, col, spotNo)   ||  
              downIs(row, col, spotNo) ||
              leftIs(row, col, spotNo) ||
              rightIs(row, col, spotNo)))  {
          return true;
        }
       if ((scratch.img[row*xs+col] == spotNo) && 
             (upIs(row, col, 0)    ||  
              downIs(row, col, 0)  ||
              leftIs(row, col, 0)  ||
              rightIs(row, col, 0)))  {
          return true;
        }
        return false;
    }

    Point[] calculateOutlineArray(int spotNo) {
       Point xy1, xy2;
       Point[] edge;
       SpotInfo s = spots[spotNo];
       if (s.x1 > 0) { xy1.x = s.x1 - 1; }
       if (s.y1 > 0) { xy1.y = s.y1 - 1; }
       if (s.x2 < (xs-1)) { 
          xy2.x = s.x2 + 1; 
       } else {
           xy2.x = s.x2;
       }
       if (s.y2 < (ys-1)) {
           xy2.y = s.y2 + 1;
       } else {
           xy2.y = s.y2;
       }
       // bounding box set, now loop through to count edge pixels
       int edgePixelCount = 0;
       int row, col;
       for (row=xy1.y; row<=xy2.y; row++) {
          for (col=xy1.x; col<=xy2.x; col++) {
             if (isEdge(row,col, spotNo)) {
                edgePixelCount++;
             }
          }
       }
       // now allot edge array, loop through to fill it with pixel locations
       edge  = new Point[edgePixelCount];
       int n = 0;
       for (row=xy1.y; row<=xy2.y; row++) {
          for (col=xy1.x; col<=xy2.x; col++) {
             if (isEdge(row, col, spotNo)) {
                edge[n++] = Point(col, row);
             }
          }
       }
       // return list of edge pixels for curve fit
       return edge;
    }

}

